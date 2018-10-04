Contamination control from shallow FASTQ
========================================

Our idea is to understand if we can perform a contamination pre-flight check from shallow data at 3x-7x coverage. We will have shallow GiaB samples data with an unknown contamination, and we can compare it to a full GiaB data known in advance. There are 2 basic approaches we want to explore: 

1. Call variants and compare AFs against a truth. Contaminated data should result in a consistently lower AF of each score, e.g. for 10% contamination, AFs should be on avarage 10% lower.
2. Count k-mers in FASTQ, the k-mer identity should be 10% less than expected for a clean sample.

- [Contamination control from shallow FASTQ](#contamination-control-from-shallow-fastq)
  - [NGSCheckMate](#ngscheckmate)
    - [FASTQ approach on clean samples](#fastq-approach-on-clean-samples)
    - [With contaminated samples](#with-contaminated-samples)
    - [Exploring VAF files](#exploring-vaf-files)
    - [Generate VAF from full BAMs](#generate-vaf-from-full-bams)
    - [Rebuilding SNPs from scratch](#rebuilding-snps-from-scratch)
    - [Rebuilding with gnomad](#rebuilding-with-gnomad)
    - [Aligning shallow BAMs](#aligning-shallow-bams)
    - [Working with new SNP list](#working-with-new-snp-list)
      - [Find VAFs from contaminated and clean FASTQ](#find-vafs-from-contaminated-and-clean-fastq)
      - [Extract VAFs from shallow BAMs](#extract-vafs-from-shallow-bams)
      - [Extract VAFs from full BAMs](#extract-vafs-from-full-bams)
      - [Compare differences based on those VAFs](#compare-differences-based-on-those-vafs)
  - [Building gnomad homozygous SNP set](#building-gnomad-homozygous-snp-set)
  - [Substituting germline variants](#substituting-germline-variants)
  - [GiaB homozygous SNPs](#giab-homozygous-snps)
  - [k-mer based approaches](#k-mer-based-approaches)
    - [-m 2 to discard single-copy](#m-2-to-discard-single-copy)
    - [-g 3234830K](#g-3234830k)
    - [Target the coverage with -c](#target-the-coverage-with--c)
    - [On contaminated data](#on-contaminated-data)
  - [ContaminationDetection tool](#contaminationdetection-tool)
  - [Conpair](#conpair)

## NGSCheckMate

To explore the first approach, we will work with [NGSCheckMate](https://github.com/parklab/NGSCheckMate) that can extract AFs from VCF, BAM, and raw FASTQ data as well. 

The tool is installed here: `cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting ; source load.sh` and requires python2.


### FASTQ approach on clean samples

Working with FASTQ directly would be useful if it allows to save time without a accuracy drop, we will start exploring this approach. NGSCheckMate uses a predefined set of ~11k SNPs which we will use to start with.

Using the following 3 shallow samples: 2 NA24631 and 1 NA24385, supposedly without contamination. Around 7x coverage each.

```
/data/cephfs/punim0010/data/FASTQ/170918_A00130_0023_BH3M2HDMXX/NA24631-PosCntl_S6_R1_001.fastq.gz  # 77,725,046 pairs ~ 3.8x * 2
/data/cephfs/punim0010/data/FASTQ/170918_A00130_0023_BH3M2HDMXX/NA24631-PosCntl_S6_R2_001.fastq.gz
/data/cephfs/punim0010/data/FASTQ/171016_A00130_0028_BH533CDMXX/NA24631_S9_R1_001.fastq.gz      	# 74,111,074 pairs ~ 3,6x * 2
/data/cephfs/punim0010/data/FASTQ/171016_A00130_0028_BH533CDMXX/NA24631_S9_R2_001.fastq.gz
/data/cephfs/punim0010/data/FASTQ/171103_A00130_0031_BH3JJYDMXX/PTC_NA24385_S1_R1_001.fastq.gz  	# 31,983,302 pairs ~ 1,6x * 2
/data/cephfs/punim0010/data/FASTQ/171103_A00130_0031_BH3JJYDMXX/PTC_NA24385_S1_R2_001.fastq.gz
```

Symlinking each into `fastqs`, creating a `fastqs.list` as per the documentation:

```
fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz       fastqs/NA24631-PosCntl_S6_R2_001.fastq.gz
fastqs/NA24631_S9_R1_001.fastq.gz       fastqs/NA24631_S9_R2_001.fastq.gz
fastqs/PTC_NA24385_S1_R1_001.fastq.gz   fastqs/PTC_NA24385_S1_R2_001.fastq.gz
```

```
python2 NGSCheckMate/ncm_fastq.py -l fastqs.list -pt NGSCheckMate/SNP/SNP.pt -O output -p 9
```

It took 47 minutes total. However internally it called `ngscheckmate_fastq` sequencially, meaning that running each sample in parallel would probably speed up the analysis:

```
/data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -p 9 -1 fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz NGSCheckMate/SNP/SNP.pt > ncm_files/NA24631-PosCntl_S6_R2_001.fastq.gz.ncm
/data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -p 9 -1 fastqs/NA24631_S9_R1_001.fastq.gz NGSCheckMate/SNP/SNP.pt > ncm_files/NA24631_S9_R2_001.fastq.gz.ncm
/data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -p 9 -1 fastqs/PTC_NA24385_S1_R1_001.fastq.gz NGSCheckMate/SNP/SNP.pt > ncm_files/PTC_NA24385_S1_R2_001.fastq.gz.ncm
Generate Data Set from output
using this bed file : NGSCheckMate/SNP/SNP.pt
```

Also interesting that it ignored the second mate in pair. For the future, we can ignore the wrapper `ncm_fastq.py` script and run `ngscheckmate_fastq` directly to better control the workflow.

The output looks correct:

```
NA24631-PosCntl_S6_R2_001.fastq.gz  matched    NA24631_S9_R2_001.fastq.gz      0.8014  1.42
NA24631-PosCntl_S6_R2_001.fastq.gz  unmatched  PTC_NA24385_S1_R2_001.fastq.gz  0.1137  0.67
NA24631_S9_R2_001.fastq.gz          unmatched  PTC_NA24385_S1_R2_001.fastq.gz  0.1174  0.67
```

We want to explore if we can downsample the shallow fastqs even further. Is 1x coverage (`-d 1`) enough?

```
time /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -1 fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz -d 1 NGSCheckMate/SNP/SNP.pt -p 9 > ncm_files/NA24631-PosCntl_S6_R2_001.D1.fastq.gz.ncm
pe=0
started reading patternfile
finished reading patternfile
finished building count array
Number of reads in the file : 77725046
Calculating subsampling rate..
read length = 149, nReads= 77725046, reference_length= 3.000000e+09, estimated depth without subsampling =3.860344, subsampling rate=0.259044, reads sampled (approximate) =20134228
Scanning fastq file..
finished reading fastq file
finished deleting hash
finished printing count array
1475.93s user 13.83s system 342% cpu 7:14.55 total
```

It took 7 minutes to run in 9 threads (however utilising only 342% cpu). Now compare with `var_ncm.py`:

```
python2 NGSCheckMate/vaf_ncm.py -f -I ncm_files -O output_clean_1sample_1x
# viewing output_clean_1sample_1x/output_all.txt
                                                                               Corr    Depth
NA24631-PosCntl_D1_R2_001.fastq.gz  matched    NA24631_S9_R2_001.fastq.gz      0.7441  0.43   # A - 1x   vs B - 3.6x
NA24631-PosCntl_D1_R2_001.fastq.gz  unmatched  PTC_NA24385_S1_R2_001.fastq.gz  0.0903  0.43   # A - 1x   vs C - 1.6x
NA24631_S9_R2_001.fastq.gz          unmatched  PTC_NA24385_S1_R2_001.fastq.gz  0.1174  0.67   # B - 3.6x vs C - 1.6x
```

The result still looks very good.

Trying now with a pre-calculated subsampling rate `0.259044` (counted in the previous `-d 1` run):

```
time /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -1 fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz --ss 0.259044 NGSCheckMate/SNP/SNP.pt -p 9 > ncm_files/NA24631-PosCntl_S6_SS025_R2_001.fastq.gz.ncm
Number of reads in the file : 77725046
1484.46s user 12.70s system 421% cpu 5:55.58 total
```

Took a bit faster, now utilising 421% cpu.

Subsampling other samples too and rerunning comparison:

```
time /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -1 fastqs/NA24631_S9_R1_001.fastq.gz --ss 0.259044 NGSCheckMate/SNP/SNP.pt -p 9 > ncm_files/NA24631_S9_SS025_R2_001.fastq.gz.ncm
Number of reads in the file : 74111074
1430.79s user 13.00s system 405% cpu 5:55.72 total

time /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -1 fastqs/PTC_NA24385_S1_R1_001.fastq.gz --ss 0.259044 NGSCheckMate/SNP/SNP.pt -p 9 > ncm_files/PTC_NA24385_S1_SS025_R2_001.fastq.gz.ncm
Number of reads in the file : 31983302
616.97s user 5.43s system 405% cpu 2:33.35 total

python2 NGSCheckMate/vaf_ncm.py -f -I ncm_files -O output_clean_1x
# should match because sample fastq:
NA24631-PosCntl_R1.fastq.gz               matched    NA24631-PosCntl_D1_R2_001.fastq.gz        0.9012  0.43   # A - 3.8x vs A - 1x
NA24631-PosCntl_D1_R2_001.fastq.gz        matched    NA24631-PosCntl_S6_SS025_R2_001.fastq.gz  0.8018  0.43   # A - 1x   vs A - 1x (re-subsampling)
NA24631_S9_R2_001.fastq.gz                matched    NA24631_S9_SS025_R2_001.fastq.gz          0.8918  0.42   # B - 3.6x vs B - 1x
PTC_NA24385_S1_R2_001.fastq.gz            matched    PTC_NA24385_S1_SS025_R2_001.fastq.gz      0.9043  0.2    # C - 1.6x vs C - 1x
# should match because same sample:
NA24631-PosCntl_R1.fastq.gz               matched    NA24631_S9_R2_001.fastq.gz                0.8014  1.42   # A - 3.8x vs B - 3.6x
NA24631-PosCntl_R1.fastq.gz               matched    NA24631_S9_SS025_R2_001.fastq.gz          0.7224  0.42   # A - 3.8x vs B - 1x
NA24631-PosCntl_D1_R2_001.fastq.gz        matched    NA24631_S9_R2_001.fastq.gz                0.7441  0.43   # A - 1x   vs B - 3.6x
NA24631-PosCntl_D1_R2_001.fastq.gz        matched    NA24631_S9_SS025_R2_001.fastq.gz          0.6707  0.42   # A - 1x   vs B - 1x
# A or B vs C - all shoud unmatch:
NA24631-PosCntl_R1.fastq.gz               unmatched  PTC_NA24385_S1_R2_001.fastq.gz            0.1137  0.67   # A - 3.8x vs C - 1.6x
NA24631-PosCntl_R1.fastq.gz               unmatched  PTC_NA24385_S1_SS025_R2_001.fastq.gz      0.1143  0.2    # A - 3.8x vs C - 1x
NA24631-PosCntl_D1_R2_001.fastq.gz        unmatched  PTC_NA24385_S1_R2_001.fastq.gz            0.0903  0.43   # A - 1x   vs C - 1.6x
NA24631-PosCntl_D1_R2_001.fastq.gz        unmatched  PTC_NA24385_S1_SS025_R2_001.fastq.gz      0.0851  0.2    # A - 1x   vs C - 1x
NA24631_S9_R2_001.fastq.gz                unmatched  PTC_NA24385_S1_R2_001.fastq.gz            0.1174  0.67   # B - 3.6x vs C - 1.6x
NA24631_S9_R2_001.fastq.gz                unmatched  PTC_NA24385_S1_SS025_R2_001.fastq.gz      0.103   0.2    # B - 3.6x vs C - 1x
NA24631_S9_SS025_R2_001.fastq.gz          unmatched  PTC_NA24385_S1_R2_001.fastq.gz            0.1069  0.42   # B - 1x   vs C - 1.6x
NA24631_S9_SS025_R2_001.fastq.gz          unmatched  PTC_NA24385_S1_SS025_R2_001.fastq.gz      0.0951  0.2    # B - 1x   vs C - 1x
```

For all comparisons, evidence is pretty strong.


### With contaminated samples

Now we want to see if we are able to detect contaminations. Out AF threshold for somatic calling is 10x, so we should be handling contamination below 10% in variant calling. So here we should try to detect contamination of 10% or above. As a test data, we emulate contaminated samples by concatenating fastq in proportion 9 to 1.

`NA24631_S9_R1_001.fastq.gz` has 74,111,074 pairs total, `PTC_NA24385_S11_R1_001.fastq.gz` has 73,623,833 pairs. We will take 90% of R1 of the former (66,699,966 reads), and 7,411,107 of the latter (to get exactly 10% of the former). This will give us the coverage of 66699966 / (3234830000 / 150) ~ 3x plus 7411107 / (3234830000 / 150) = 0.34x. 

We assume that the reads in fastq are placed randomly (which should be safe for us), so can use `head` for subsampling:

```
time gunzip -c NA24631_S9_R1_001.fastq.gz      | head -n $(py '66699966 * 4') | gzip -c > NA24631_S9_R1_001.head_3x.fastq.gz
time gunzip -c PTC_NA24385_S11_R1_001.fastq.gz | head -n $(py '7411107 * 4')  | gzip -c > PTC_NA24385_S11_R1_001.head_034x.fastq.gz

# gunzip -c NA24631_S9_R2_001.fastq.gz          184.31s user 11.49s system  7% cpu 44:35.48 total
# head -n $(py '66699966 * 4')                   25.76s user 27.94s system  2% cpu 44:26.05 total
# gzip -c >> NA24631_S9.head_6x.fastq.gz       2580.25s user  8.86s system 97% cpu 44:26.06 total
```

The coverage might be quite low, so we will also test 6x+0.7x (129393200 + 15095873) by using the R2 reads as well:

```
# 6x + 0.7x:
cp             NA24631_S9_R1_001.head_3x.fastq.gz                                          NA24631_S9.head_6x.fastq.gz
time gunzip -c NA24631_S9_R2_001.fastq.gz      | head -n $(py '66699966 * 4') | gzip -c >> NA24631_S9.head_6x.fastq.gz
cp             PTC_NA24385_S11_R1_001.head_034x.fastq.gz                                   PTC_NA24385_S11.head_07x.fastq.gz
time gunzip -c PTC_NA24385_S11_R2_001.fastq.gz | head -n $(py '7411107 * 4')  | gzip -c >> PTC_NA24385_S11.head_07x.fastq.gz

# gunzip -c PTC_NA24385_S11_R2_001.fastq.gz      21.28s user  1.31s system  6% cpu  5:29.70 total
# head -n $(py '7411107 * 4')                     2.81s user  3.03s system  1% cpu  5:20.03 total
# gzip -c >> PTC_NA24385_S11.head_07x.fastq.gz  311.62s user  1.11s system 97% cpu  5:20.04 total
```

For the sake of experiment, aslo trying to subsample correctly just in case:

```
time seqtk sample -s 11 PTC_NA24385_S11_R1_001.fastq.gz 0.9 | gzip -c > PTC_NA24385_S11_R1_001.seqtk_90pct.fastq.gz
# seqtk:    343.89s user 20.75s system 10% cpu 1:00:30.13 total
# gzip:     3100.37s user 13.22s system 85% cpu 1:00:30.15 total
# Timing for seqtk only:    330.30s user 17.57s system 61% cpu 9:24.57 total
```

It seem to be taking quite a lot of CPU time, though, and shouldn't affect the results.

Merging:

```
cat NA24631_S9_R1_001.head_3x.fastq.gz PTC_NA24385_S11_R1_001.head_034x.fastq.gz > NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz
cat NA24631_S9.head_6x.fastq.gz        PTC_NA24385_S11.head_07x.fastq.gz         > NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz
```

Running NGSCheckMate:

```
time /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -p 9 -1 fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz NGSCheckMate/SNP/SNP.pt > ncm_files/NA24631_S9__3x__PTC_NA24385_S11__034x.ncm
# Number of reads in the file : 74 111 073
# 3690.93s user 7.44s system 640% cpu 9:37.36 total
 
time /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -p 9 -1 fastqs/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz NGSCheckMate/SNP/SNP.pt > ncm_files/NA24631_S9__6x__PTC_NA24385_S11__07x.ncm
# Number of reads in the file : 150 074 923
# 7287.73s user 17.90s system 650% cpu 18:42.33 total

python2 NGSCheckMate/vaf_ncm.py -f -I ncm_files -O output_contaminated

# vs. clean B
NA24631-PosCntl_R1.fastq.gz            matched    NA24631_S9_R2_001.fastq.gz             0.8014    # A - 3.8x  vs.  B - 3.6x
NA24631-PosCntl_D1_R2_001.fastq.gz     matched    NA24631_S9_R2_001.fastq.gz             0.7441    # A - 1x    vs.  B - 3.6x 
PTC_NA24385_S1_R2_001.fastq.gz         unmatched  NA24631_S9_R2_001.fastq.gz             0.1174    # C - 1.6x  vs.  B - 3.6x
PTC_NA24385_S1_SS025_R2_001.fastq.gz   unmatched  NA24631_S9_R2_001.fastq.gz             0.103     # C - 1x    vs.  B - 3.6x

# vs. 3x+0.34x contaminated
NA24631-PosCntl_D1_R2_001.fastq.gz     matched    NA24631_S9__3x__PTC_NA24385_S11__034x  0.7113    # A - 3.8x  vs.  B - 3x + C - 0.34x
NA24631-PosCntl_R1.fastq.gz            matched    NA24631_S9__3x__PTC_NA24385_S11__034x  0.7668    # A - 1x    vs.  B - 3x + C - 0.34x
NA24631_S9_R2_001.fastq.gz             matched    NA24631_S9__3x__PTC_NA24385_S11__034x  0.9546    # B - 3.6x  vs.  B - 3x + C - 0.34x
NA24631_S9_SS025_R2_001.fastq.gz       matched    NA24631_S9__3x__PTC_NA24385_S11__034x  0.8621    # B - 1x    vs.  B - 3x + C - 0.34x
PTC_NA24385_S1_R2_001.fastq.gz         unmatched  NA24631_S9__3x__PTC_NA24385_S11__034x  0.1837    # C - 1.6x  vs.  B - 3x + C - 0.34x  
PTC_NA24385_S1_SS025_R2_001.fastq.gz   unmatched  NA24631_S9__3x__PTC_NA24385_S11__034x  0.1719    # C - 1x    vs.  B - 3x + C - 0.34x

# vs. 6x+0.7x contaminated
NA24631-PosCntl_D1_R2_001.fastq.gz     matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.6224    # A - 3.8x  vs.  B - 6x + C - 0.7x
NA24631-PosCntl_R1.fastq.gz            matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.6793    # A - 1x    vs.  B - 6x + C - 0.7x
NA24631_S9_R2_001.fastq.gz             matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.816     # B - 3.6x  vs.  B - 6x + C - 0.7x
NA24631_S9_SS025_R2_001.fastq.gz       matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.7357    # B - 1x    vs.  B - 6x + C - 0.7x
PTC_NA24385_S1_R2_001.fastq.gz         unmatched  NA24631_S9__6x__PTC_NA24385_S11__07x   0.1488    # C - 1.6x  vs.  B - 6x + C - 0.7x
PTC_NA24385_S1_SS025_R2_001.fastq.gz   unmatched  NA24631_S9__6x__PTC_NA24385_S11__07x   0.1403    # C - 1x    vs.  B - 6x + C - 0.7x

# contaminated vs contaminated
NA24631_S9__3x__PTC_NA24385_S11__034x  matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.9359    # B - 3x + C - 0.34x  vs B - 6x + C - 0.7x
```

We see that 6x+0.7x yields a higher similarity with PTC_NA24385 (~0.14) than when comparing against a pure NA24631 (~0.11), and at the same time it yields a lower similarity with NA24631 (~0.6) than when comparing against a pure NA24631 (~0.8) - which looks promising, however doesn't give us exact thresholds to confidently tell one from the other. We will dig into the VAF files to see if we can extract the contamination rate from them - e.g., calculate a median differece between AFs for all matching calls.


### Exploring VAF files

We calculate the difference between AFs. We exclude lines with `NA`.

```
# Identical samples:
paste NA24631_S9_R2_001.fastq.gz.ncm NA24631-PosCntl_S6_R2_001.fastq.gz.ncm    | grep -v NA | awk 'BEGIN {OFS="\t"} function abs(v) {return v < 0 ? -v : v} { sum += abs($4 - $8) } END { print sum / NR }'
# 0.13942

# Different samples:
paste NA24631_S9_R2_001.fastq.gz.ncm PTC_NA24385_S11_R1_001.fastq.gz.ncm       | grep -v NA | awk 'BEGIN {OFS="\t"} function abs(v) {return v < 0 ? -v : v} { sum += abs($4 - $8) } END { print sum / NR }'
# 0.444964

# Contaminated - 3.5x:
paste NA24631_S9_R2_001.fastq.gz.ncm NA24631_S9__3x__PTC_NA24385_S11__034x.ncm | grep -v NA | awk 'BEGIN {OFS="\t"} function abs(v) {return v < 0 ? -v : v} { sum += abs($4 - $8) } END { print sum / NR }'
# 0.0552495

# Contaminated - 7x:
paste NA24631_S9_R2_001.fastq.gz.ncm NA24631_S9__6x__PTC_NA24385_S11__07x.ncm  | grep -v NA | awk 'BEGIN {OFS="\t"} function abs(v) {return v < 0 ? -v : v} { sum += abs($4 - $8) } END { print sum / NR }'
# 0.0899408
```

Obviously, it's not sensitive enough to detect contamintaion. However, in reality we will be comparing the shallow QC datasets against full clean pre-exusting data.


### Generate VAF from full BAMs

We will generate a VAF file from a full BAM and explore how our shallow VAF files compare to it.

NGSCheckMate can run from VCFs, so trying to generate a VAF from VCF:

```
bcftools view /data/cephfs/punim0010/projects/Hsu_WGS_Validation/WGS-GiaB-merged/final/2017-11-13_giab-merged/HanChinese-1-ensemble-annotated.vcf.gz -R NGSCheckMate/SNP/SNP_GRCh37_hg19_woChr.bed > NA24631_full.vcf
bgzip NA24631_full.vcf
tabix -p vcf NA24631_full.vcf.gz

# convert VCF to NA24631_full.ncm VCF file. For each record in SNP.bed, tabix one from VCF. Or just vardict using BED file?
# bcftools query -f "X\t%REF\t%ALT\t%AF\n" NA24631_full.vcf | tsv
index   ref     alt     vaf
0       0       0       NA
1       NA      NA      NA
2       0       2       1.000000

samtools mpileup -auf /data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -l NGSCheckMate/SNP/SNP_GRCh37_hg19_woChr.bed /data/cephfs/punim0010/projects/Hsu_WGS_Validation/WGS-GiaB-merged/final/NA24631-1KC/NA24631-1KC-ready.bam -t DP,AD | bcftools view -H > NA24631_mpileup_a.vcf
```

Internally, for VCF input ncm.py runs same samtools mpileup - so we can do it manually as well.

```
python2 NGSCheckMate/ncm.py -B -l vcfs.list -O ncm_from_vcfs -N ncm_from_vcfs -bed NGSCheckMate/SNP/SNP_GRCh37_hg19_woChr.bed
# it runs internally the following:
# samtools mpileup -uf /data/cephfs/punim0010/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -l NGSCheckMate/SNP/SNP_GRCh37_hg19_woChr.bed /data/cephfs/punim0010/projects/Hsu_WGS_Validation/WGS-GiaB-merged/final/2017-11-13_giab-merged/AshkenazimJew-2-ensemble-annotated.vcf.gz | bcftools call -c > ncm_from_vcfs/AshkenazimJew-2-ensemble-annotated.vcf
```

However, the SNPs BED file looks fishy - some SNPs are 45 bases long https://github.com/parklab/NGSCheckMate/issues/15

```
cat NGSCheckMate/SNP/SNP_GRCh37_hg19_wChr.bed | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $3-$2 }' | sort -k6,6 | less
```

All those rsIDs are deleted from the NCBI database, as well as many others in the list:

```
vcfanno dbsnp.toml ncm_from_bams/NA24631-1KC-ready.vcf > ncm_from_bams/NA24631-1KC-ready.anno.vcf
grep rs ncm_from_bams/NA24631-1KC-ready.anno.vcf | wc
  10934  109323 2258957
grep -v rs ncm_from_bams/NA24631-1KC-ready.anno.vcf | wc
  10459  103828 1451221
```

So we would need to rebuild the SNP list.


### Rebuilding SNPs from scratch

We will re-build the SNPs from scratch. We will prepare a list of SNPs with MAF close to 50%, and run NGSCheckMate's [patterngenerator](https://github.com/parklab/NGSCheckMate#1-patterngenerator). "It requires a bed file containing a set of SNP positions, a genome reference file (both FASTA and bowtie 1 index) and the bowtie alignment program (http://bowtie-bio.sourceforge.net/index.shtml)." In order to make a BED file, "We required uniqueness of the k-mer sequence flanking the SNPs and used 11 696 SNPs (more details below). Our simulation showed that there was no difference in the distribution of VAF correlations between using the 21 067 and the reduced 11 696 SNP sets."

We will overlap with a list of single-allelic dbSNP SNPs in confident regions with GMAF>10%, and re-build the pt file from that:

```
cd new_SNPs
wget https://github.com/AstraZeneca-NGS/ClearUp/blob/master/clearup/snps/dbsnp_maf10pct.hg19.no_selfchain_gc25-30_65-70_lowcomp50.autosomal.bed.gz\?raw\=true -O dbsnp_maf10pct.hg19.no_selfchain_gc25-30_65-70_lowcomp50.autosomal.bed.gz

bedtools sort -i ../NGSCheckMate/SNP/SNP_GRCh37_hg19_wChr.bed > SNP_hg19.sorted.bed

bedtools intersect -a SNP_hg19.sorted.bed.gz -b dbsnp_maf10pct.hg19.no_selfchain_gc25-30_65-70_lowcomp50.autosomal.bed.gz > SNP_hg19.sorted.maf10.bed

cat SNP_hg19.sorted.maf10.bed | py -x "x[3:]" > SNP_GRCh37.sorted.maf10.bed

perl ../NGSCheckMate/patterngenerator/makesnvpattern.pl SNP_GRCh37.sorted.maf10.bed /data/cephfs/punim0010/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa /data/cephfs/punim0010/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/bowtie/GRCh37.fa . new_SNPs
```

### Rebuilding with gnomad

```

cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/gnomad_SNPs
bcftools query gnomad.genomes.r2.0.2.sites.coding_only.chr1-22.SNPS.SINGLE.PASS.COMMON.vcf.bgz -f "%CHROM\t%POS\t%POS\t%ID\n" | awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4}' > gnomad.genomes.r2.0.2.sites.coding_only.chr1-22.SNPS.SINGLE.PASS.COMMON.bed
../NGSCheckMate/patterngenerator/makesnvpattern.pl gnomad.genomes.r2.0.2.sites.coding_only.chr1-22.SNPS.SINGLE.PASS.COMMON.bed /data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa /data/cephfs/punim0010/extras/vlad/bcbio/genomes/Hsapiens/GRCh37/bowtie/GRCh37.fa . gnomad_SNP
```


### Aligning shallow BAMs

We also want to see how sensitive the FASTQ method is compared to calling variants from aligned BAMs. We running bcbio for positive controls:

```
# Preparing left reads in addition to NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz:
ln -s NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz NA24631_S9__3x__PTC_NA24385_S11__034x_1.fastq.gz
time gunzip -c NA24631_S9_R2_001.fastq.gz      | head -n $(py '66699966 * 4') | gzip -c > NA24631_S9_R2_001.head_3x.fastq.gz
time gunzip -c PTC_NA24385_S11_R2_001.fastq.gz | head -n $(py '7411107 * 4')  | gzip -c > PTC_NA24385_S11_R2_001.head_034x.fastq.gz
cat NA24631_S9_R2_001.head_3x.fastq.gz PTC_NA24385_S11_R2_001.head_034x.fastq.gz > NA24631_S9__3x__PTC_NA24385_S11__034x_2.fastq.gz

cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bcbio

bcbio_nextgen.py -w template std_workflow_positiveControl.yaml bcbio.csv /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x_* /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/fastqs/NA24631_S9_R?_001.fastq.gz
```

Compare with FASTQ methods:

```
paste NA24631_S9.ncm ../ncm_files/NA24631_S9_R2_001.fastq.gz.ncm | grep -v NA | awk 'BEGIN {OFS="\t"; print "", "", "", "BAM", "", "", "", "FASTQ", "BAM - FASTQ"} function abs(v) {return v < 0 ? -v : v} { print $1, $2, $3, $4, $5, $6, $7, $8, abs($4 - $8) }' | tsv
```

### Working with new SNP list

#### Find VAFs from contaminated and clean FASTQ

```
NCM=/data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq
PT=new_SNPs/new_SNPs.pt
time $NCM -p 10 -1 fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz $PT > ncm_files/NA24631_S9__3x__PTC_NA24385_S11__034x.ncm &
time $NCM -p 10 -1 fastqs/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz  $PT > ncm_files/NA24631_S9__6x__PTC_NA24385_S11__07x.ncm &
time $NCM -p 10 -1 fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz             $PT > ncm_files/NA24631-PosCntl_S6_R2_001.fastq.gz.ncm &
time $NCM -p 10 -1 fastqs/NA24631_S9_R1_001.fastq.gz                     $PT > ncm_files/NA24631_S9_R2_001.fastq.gz.ncm &
time $NCM -p 10 -1 fastqs/PTC_NA24385_S11_R1_001.fastq.gz                $PT > ncm_files/PTC_NA24385_S11_R1_001.fastq.gz.ncm &

python2 NGSCheckMate/vaf_ncm.py -f -I ncm_files -O output_contaminated

```

#### Extract VAFs from shallow BAMs

```
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bcbio/bcbio/final/NA24631_S9__3x__PTC_NA24385_S11__034x
mkdir mpileup_work
REF=/data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
SNPS=/data/cephfs/punim0010/projects/Saveliev_Fingerprinting/new_SNPs/SNP_GRCh37.sorted.maf10.bed
cat $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs -p 31 "samtools mpileup -uf $REF -r {} NA24631_S9__3x__PTC_NA24385_S11__034x-ready.bam -a -t DP,AD | bcftools view -H -o mpileup_work/NA24631_S9__3x__PTC_NA24385_S11__034x-ready.bam.mpileup_{}"

ls -1 mpileup_work/NA24631_S9__3x__PTC_NA24385_S11__034x-ready.bam.mpileup_* | parallel -k "if [ -s {} ] ; then cat {} ; else echo "NA {}" ; fi" > NA24631_S9__3x__PTC_NA24385_S11__034x-mpileup_all

sort --general-numeric-sort -k1,1 -k2,2 NA24631_S9__3x__PTC_NA24385_S11__034x-mpileup_all | uniq > NA24631_S9__3x__PTC_NA24385_S11__034x-mpileup_all_sorted

# In result, getting lines as follows - NA for uncalled and proper VCF records for called:
# 1       168013850       .       T       <*>     0       .       DP=2;I16=2,0,0,0,74,2738,0,0,120,7200,0,0,35,725,0,0;QS=1,0;MQ0F=0      PL:DP:ADF:AD    0,6,67:2:2,0:2,0
# NA NA24631_S9__3x__PTC_NA24385_S11__034x-sort.mpileup_10:100018844-100018844
```

Now for all records, generate ncm file like

```
index   ref     alt     vaf
0       0       0       NA
1       NA      NA      NA
2       0       2       1.000000
```

Add index, and if it's NA - add NA, otherwise add `AD[1] / sum(AD)`

```
python make_vaf.py NA24631_S9__3x__PTC_NA24385_S11__034x-mpileup_all_sorted > NA24631_S9__3x__PTC_NA24385_S11__034x.ncm
```

Same for non-contaminated file:

```
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bcbio/bcbio/final/NA24631_S9
mkdir mpileup_work
cat $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs -p 31 "samtools mpileup -uf $REF -r {} NA24631_S9-ready.bam -a -t DP,AD | bcftools view -H -o mpileup_work/NA24631_S9-ready.bam.mpileup_{}"
ls -1 mpileup_work/NA24631_S9-ready.bam.mpileup_* | parallel -k "if [ -s {} ] ; then cat {} ; else echo "NA {}" ; fi" > NA24631_S9-mpileup_all
sort --general-numeric-sort -k1,1 -k2,2 NA24631_S9-mpileup_all | uniq > NA24631_S9-mpileup_all_sorted
python make_vaf.py NA24631_S9-mpileup_all_sorted > NA24631_S9.ncm
```

And the contamination sample:

```
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bcbio/bcbio/final/PTC_NA24385_S11
mkdir mpileup_work
cat $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs -p 31 "samtools mpileup -uf $REF -r {} PTC_NA24385_S11-ready.bam -a -t DP,AD | bcftools view -H -o mpileup_work/PTC_NA24385_S11-ready.bam.mpileup_{}"
ls -1 mpileup_work/PTC_NA24385_S11-ready.bam.mpileup_* | parallel -k "if [ -s {} ] ; then cat {} ; else echo "NA {}" ; fi" > PTC_NA24385_S11-mpileup_all
sort --general-numeric-sort -k1,1 -k2,2 PTC_NA24385_S11-mpileup_all | uniq > PTC_NA24385_S11-mpileup_all_sorted
python make_vaf.py PTC_NA24385_S11-mpileup_all_sorted > PTC_NA24385_S11.ncm
```


#### Extract VAFs from full BAMs

```
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/ncm_from_bams/mpileups
cat $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs -p 31 "samtools mpileup -uf $REF -r {} /data/cephfs/punim0010/projects/Hsu_WGS_Validation/WGS-GiaB-merged/final/NA24631-1KC/NA24631-1KC-ready.bam -a -t DP,AD | bcftools view -H -o NA24631-1KC-ready.mpileup/NA24631-1KC.mpileup_{}"
ls -1 NA24631-1KC-ready.mpileup/NA24631-1KC.mpileup_* | parallel -k "if [ -s {} ] ; then cat {} ; else echo "NA {}" ; fi" > NA24631-1KC-ready.mpileup/NA24631-1KC.mpileup_all
sort --general-numeric-sort -k1,1 -k2,2 NA24631-1KC-ready.mpileup/NA24631-1KC.mpileup_all | uniq > NA24631-1KC.mpileup_all_sorted
cd ..
python make_vaf.py mpileups/NA24631-1KC.mpileup_all_sorted > NA24631-1KC.ncm
```

#### Compare differences based on those VAFs

```
# Full - shallow-clean:
paste NA24631-1KC.ncm NA24631_S9.ncm                            | grep -v NA | awk 'BEGIN {OFS="\t"} function abs(v) {return v < 0 ? -v : v} { sum += abs($4 - $8) } END { print sum / NR }'
# 0.0751213

# Full - shallow-contaminated:
paste NA24631-1KC.ncm NA24631_S9__3x__PTC_NA24385_S11__034x.ncm | grep -v NA | awk 'BEGIN {OFS="\t"} function abs(v) {return v < 0 ? -v : v} { sum += abs($4 - $8) } END { print sum / NR }'
# 0.102715

# Full - shallow-contamination-alone:
paste NA24631-1KC.ncm PTC_NA24385_S11.ncm                       | grep -v NA | awk 'BEGIN {OFS="\t"} function abs(v) {return v < 0 ? -v : v} { sum += abs($4 - $8) } END { print sum / NR }'
# 0.415981
```

To explore in detail, use:

```
paste NA24631-1KC.mpileup_all.ncm NA24631_S9.ncm NA24631_S9__3x__PTC_NA24385_S11__034x.ncm PTC_NA24385_S11.ncm | grep -v NA | grep -v vaf | awk 'BEGIN {OFS="\t"; print "A", "Shal A", "Shal A w/B", "Shal B", "A - Shal A", "A - Shal A w/B", "A - Shal B"} function abs(v) {return v < 0 ? -v : v} { print $4, $8, $12, $16, abs($4 - $8), abs($4 - $12), abs($4 - $16) }' | tsv
# returns: 
```

## Building gnomad homozygous SNP set

If we select variants that predominantly appear as homozygous in the population, evenly distributed between HOM REF and HOM ALT, then any non-homozygous call will indicate contamination. We prepare a `gnomad` file with such SNPs:

We have the bcbio `gnomad` file:

```
bcftools view -H /data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/variation/gnomad_genome.vcf.gz | wc
260,657,435
```

And we also can download the most recent one from the website:

```
http://gnomad.broadinstitute.org/downloads
# Coding only:
bcftools view -H gnomad.genomes.r2.0.2.sites.coding_only.chr1-22.vcf.bgz | wc
4,500,535
```

First, we select common (20% <= AF <= 80%) single-allelic PASSing SNPs:

```
GNO=gnomad.genomes.r2.0.2.sites.coding_only.chr1-22.vcf.bgz
GNO_FILT=gnomad.genomes.r2.0.2.sites.coding_only.chr1-22.SNPS.SINGLE.PASS.COMMON.vcf.bgz

bcftools view -f.,PASS --max-alleles 2 -v snps $GNO -Ob | bcftools filter -i "AN > 1000 && AC/AN <= 0.8 && AC/AN >= 0.2" -Oz -o $GNO_FILT
# 40481  # wc $GNO_FILT 
```

Getting 40481 variants.

Now we want to see how many there are SNPs with Hom/Het ratio above 5. Gnomad reports AN (total sample count), AC (mutant sample count), and Hom (homozygouys mutant sample count). We can calculate Het ratio as follows:
AN = AC + REF_Hom => REF_Hom = AN - AC
Total_hom = Hom + REF_Hom = Hom + AN - AC
AC = Het + Hom => Het = AC - Hom
Total_hom / Het = (Hom + AN - AC) / (AC - Hom)

```
bcftools filter -i "(Hom + AN - AC) / (AC - Hom) > 5" $GNO_FILT | bcftools query -f "%AN \t %AC \t %Hom \n" | grep -v ',' | awk 'BEGIN {OFS="\t"; print "AN", "AC", "Hom", "AF", "Het", "Total_hom", "Het / Total_hom"} { print $1, $2, $3, $2/$1, $2-$3, $3+$1-$2, ($3+$1-$2) / ($2-$3)}'
# 125  # wc
```

And we are getting only 125 variants. Also, all of those variants have only 20% of AF, which is the reason they have a low Het ratio - most of the cases are REF Hom. We want the signal to be stronger, so might require Het/Hom ratio to be calculated only from ALT:
AC = Het + Hom => Het = AC - Hom
Hom / Het = Hom / (AC - Hom)

```
bcftools filter -i "Hom / (AC - Hom) > 5" $GNO_FILT | bcftools query -f "%AN \t %AC \t %Hom \n" | grep -v ',' | awk 'BEGIN {OFS="\t"; print "Total", "Mut", "MutHom", "MutFreq", "MutHet=Mut-HomHom", "MutHet / MutHom", "TotalHom=Total-Mut+MutHom", "HomFreq=TotalHom/Total", "HetFreq"} { if ($2/$1 > 0.3 && $2/$1 < 0.7) { print $1, $2, $3, $2/$1, $2-$3, $3 / ($2-$3), $1-$2+$3, ($1-$2+$3)/$1, ($2-$3) /
$2 } }' | tsv
# 0  # wc
```

Getting 0 variants. It might come from the fact that we are using only coding SNPs. Downloading the full gnomad database and rebuilding the variant set:

```
GNO_FILT=/data/cephfs/punim0010/projects/Saveliev_Fingerprinting/gnomad/gnomad.genomes.r2.0.2.sites.snps_pass_common.vcf.bgz
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/gnomad
parallel -j 22 "echo {} && wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr{}.vcf.bgz ." ::: seq 1 22
parallel -j 22 "echo {} ; bcftools view -f.,PASS --max-alleles 2 -v snps {} -Ob | bcftools filter -i 'AN > 1000 && AC/AN <= 0.8 && AC/AN >= 0.2' -Oz -o {.}.snps_pass_common.vcf.bgz" ::: *chr*.vcf.bgz
parallel -j 22 "echo {} ; tabix -p vcf {}" ::: *.snps_pass_common.vcf.bgz
bcftools merge *.snps_pass_common.vcf.bgz -Oz -o $GNO_FILT
# 3060191  # wc $GNO_FILT 
```

Getting 3060191 common variants. Now finding predominantly homozygous variants. In chr21:

```
bcftools filter -i "Hom / (AC - Hom) > 5" gnomad.genomes.r2.0.2.sites.chr21.vcf.snps_pass_common.vcf.bgz | bcftools query -f "%AN \t %AC \t %Hom \n" | grep -v ',' | awk 'BEGIN {OFS="\t"; print "AN", "AC", "Hom", "AF", "Het", "Het / Hom"} { print $1, $2, $3, $2/$1, $2-$3, $3+$1-$2, $3 / ($2-$3)}'
```

In chr21, there are 0.



## Substituting germline variants

Now we experiment with calling germline variants in shallow BAM files.

Test:

```
grep 125479363 $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs -p 31 "samtools mpileup -uf $REF -r {} PTC_NA24385_S11-ready.bam -a -t DP,AD | bcftools call -c | bcftools view -H"
```

```
REF=/data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
SNPS=/data/cephfs/punim0010/projects/Saveliev_Fingerprinting/new_SNPs/SNP_GRCh37.sorted.maf10.bed

SAMPLE=NA24631_S9__3x__PTC_NA24385_S11__034x
BAM=$SAMPLE-ready.bam
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bcbio/bcbio/final/$SAMPLE
#cat $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs --ordered -p 31 "samtools mpileup -uf $REF -r {} $BAM -a -t DP,AD | bcftools call -c | bcftools view -Oz -o mpileup_work/$SAMPLE-ready-{}.vcf && tabix mpileup_work/$SAMPLE-ready-{}.vcf.gz"
#bcftools merge mpileup_work/$SAMPLE-ready-*.vcf.gz -o $SAMPLE-ready.vcf
samtools mpileup -uf $REF -l $SNPS $BAM -a -t DP,AD | bcftools call -c | bcftools view -Oz -o $SAMPLE-ready.vcf.gz && tabix $SAMPLE-ready.vcf.gz

SAMPLE=NA24631_S9
BAM=$SAMPLE-ready.bam
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bcbio/bcbio/final/$SAMPLE
samtools mpileup -uf $REF -l $SNPS $BAM -a -t DP,AD | bcftools call -c | bcftools view -Oz -o $SAMPLE-ready.vcf.gz && tabix $SAMPLE-ready.vcf.gz

SAMPLE=PTC_NA24385_S11
BAM=$SAMPLE-ready.bam
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bcbio/bcbio/final/$SAMPLE
samtools mpileup -uf $REF -l $SNPS $BAM -a -t DP,AD | bcftools call -c | bcftools view -Oz -o $SAMPLE-ready.vcf.gz && tabix $SAMPLE-ready.vcf.gz

cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/ncm_from_bams/mpileups
SAMPLE=NA24631-1KC
BAM=/data/cephfs/punim0010/projects/Hsu_WGS_Validation/WGS-GiaB-merged/final/NA24631-1KC/NA24631-1KC-ready.bam
~/bin/sambamba slice $BAM -L $SNPS > $SAMPLE-ready.SNPS.bam
BAM=$SAMPLE-ready.SNPS.bam
samtools mpileup -uf $REF -l $SNPS $BAM -a -t DP,AD | bcftools call -c | bcftools view -Oz -o $SAMPLE-ready.vcf.gz && tabix $SAMPLE-ready.vcf.gz
#cat $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs --ordered -p 31 "samtools mpileup -uf $REF -r {} $BAM -a -t DP,AD | bcftools call -c | bcftools view -Oz -o mpileup_work/$SAMPLE-ready-{}.vcf && tabix mpileup_work/$SAMPLE-ready-{}.vcf.gz"
```


## GiaB homozygous SNPs

We don't need gnomad SNPs if we always know the samples being studied for contamination. We can instead use the full BAM and extract clean homozugous SNPs, and check if they keep being clean homozygous in contaminated shallow data.

Extracting hom from NA24631 germline vcf:

```
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/hom
bcftools view -f.,PASS -g hom -v snps /data/cephfs/punim0010/data/Results_Local/Accreditation/GiaB/2017-11-13/WGS-GiaB-merged/final/2017-11-13_giab-merged/HanChinese-2-ensemble-annotated.vcf.gz -Oz -o NA24631-2KC.hom.snps.vcf.gz
bcftools filter -i "FORMAT/DP>=30 && FORMAT/AF==1" NA24631-1KC.hom.snps.vcf.gz -Oz -o NA24631-1KC.hom.snps.dp30af1.vcf.gz
# 1028312 variants
gunzip -c NA24631-1KC.hom.snps.dp30af1.vcf.gz | bioawk -tc vcf '{ print $chrom, $pos-1, $pos }' > NA24631-1KC.hom.snps.dp30af1.bed
```

Now call pileups:

```
mkdir mpileup_work
REF=/data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
SNPS=NA24631-1KC.hom.snps.dp30af1.bed
cat $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs -p 31 "samtools mpileup -uf $REF -r {} ../bams/NA24631_S9__3x__PTC_NA24385_S11__034x-sort.bam -a -t DP,AD | bcftools view -H -o mpileup_work/NA24631_S9__3x__PTC_NA24385_S11__034x.mpileup_{}"

for f in mpileup_work/NA24631_S9__3x__PTC_NA24385_S11__034x.mpileup_* ; do if [ -s $f ] ; then cat $f ; else echo "NA $f" ; fi ; done > NA24631_S9__3x__PTC_NA24385_S11__034x-mpileup_all

sort --general-numeric-sort -k1,1 -k2,2 NA24631_S9__3x__PTC_NA24385_S11__034x-mpileup_all | uniq > NA24631_S9__3x__PTC_NA24385_S11__034x-mpileup_all_sorted

# In result, getting lines as follows - NA for uncalled and proper VCF records for called:
# 1       168013850       .       T       <*>     0       .       DP=2;I16=2,0,0,0,74,2738,0,0,120,7200,0,0,35,725,0,0;QS=1,0;MQ0F=0      PL:DP:ADF:AD    0,6,67:2:2,0:2,0
# NA NA24631_S9__3x__PTC_NA24385_S11__034x-sort.mpileup_10:100018844-100018844
```

Now for all records, generate ncm file like:

```
index   ref     alt     vaf
0       0       0       NA
1       NA      NA      NA
2       0       2       1.000000
```

Add index, and if it's NA - add NA, otherwise add `AD[1] / sum(AD)`

```
python make_vaf.py NA24631_S9__3x__PTC_NA24385_S11__034x-mpileup_all_sorted > NA24631_S9__3x__PTC_NA24385_S11__034x.ncm
```

Same for non-contaminated file:

```
REF=/data/cephfs/punim0010/local/stable/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
SNPS=NA24631-1KC.hom.snps.dp30af1.bed
cat $SNPS | awk '{print $1":"($2+1)"-"$3}' | gargs -p 31 "samtools mpileup -uf $REF -r {} ../bams/NA24631_S9-sort.bam -a -t DP,AD | bcftools view -H -o mpileup_work/NA24631_S9.mpileup_{}"
for f in mpileup_work/NA24631_S9.mpileup_* ; do if [ -s $f ] ; then cat $f; else echo "NA $f" ; fi ; done > NA24631_S9-mpileup_all
sort --general-numeric-sort -k1,1 -k2,2 NA24631_S9-mpileup_all | uniq > NA24631_S9-mpileup_all_sorted
python make_vaf.py NA24631_S9-mpileup_all_sorted > NA24631_S9.ncm
```

We expect those varaints to be close to AF=1 on average in full sample, and below AF in contaminated sample:

```
cut -f4 NA24631_S9__3x__PTC_NA24385_S11__034x__100.ncm | awk '{ sum += $1 } END { print sum/NR }'
cut -f4 NA24631_S9.ncm | awk '{ sum += $1 } END { print sum/NR }'
```

For first 100 snps, the result is 0.961553, 0.955937 for contaminated.

For 100k snps, 0.993578 for clean one, 0.956101 for contaminated. Actually same numbers we get for 1k snps, and it's stable on any number of snps (checked up until 1m).


## k-mer based approaches

We also explored another way to match samples - a k-mer based [mash](http://mash.readthedocs.io/en/latest/) approach. Installing:

```
conda install -c bioconda mash -y
```

Running in 9 threads (make sure to use `-r` to specify that the input is reads!):

```
time mash sketch fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz -r -p 9 -o mash/NA24631-PosCntl.msh
Estimated genome size: 2.40473e+09
Estimated coverage:    3.343
954.52s user 1.84s system 99% cpu 15:56.38 total

time mash sketch fastqs/NA24631_S9_R1_001.fastq.gz -r -p 4 -o mash/NA24631.msh &
Estimated genome size: 2.51025e+09
Estimated coverage:    3.094
917.82s user 2.10s system 99% cpu 15:19.93 total

time mash sketch fastqs/PTC_NA24385_S1_R1_001.fastq.gz -r -p 4 -o mash/PTC_NA24385.msh &
Estimated genome size: 1.80646e+09
Estimated coverage:    1.864
395.67s user 0.85s system 99% cpu 6:36.56 total
```

Calculating the distance with `mash dist`:

```
fastqs/NA24631_S9			fastqs/NA24631-PosCntl   0.00932572      0       698/1000
fastqs/NA24631_S9			fastqs/PTC_NA24385       0.0161637       0       553/1000
fastqs/NA24631-PosCntl_S6	fastqs/PTC_NA24385       0.014871        0       577/1000
```

The related samples have the distance almost 700, while unrelated are below 600. However the signal doesn't look strong.

### -m 2 to discard single-copy

Now using `-m 2` to discard single-copy kmers:

```
time mash sketch fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz -r -p 9 -m 2 -o mash/NA24631-PosCntl.m2.msh
Estimated genome size: 1.81814e+09
Estimated coverage:    4.27
957.59s user 1.79s system 99% cpu 15:59.47 total

time mash sketch fastqs/NA24631_S9_R1_001.fastq.gz -r -p 9 -m 2 -o mash/NA24631.m2.msh
Estimated genome size: 1.83433e+09
Estimated coverage:    4.025
902.89s user 2.16s system 99% cpu 15:05.26 total

time mash sketch fastqs/PTC_NA24385_S1_R1_001.fastq.gz -r -p 9 -m 2 -o mash/PTC_NA24385.m2.msh
Estimated genome size: 9.03879e+08
Estimated coverage:    2.974
377.47s user 1.66s system 99% cpu 6:21.62 total
```

The difference is more clear - unrelated samples are below 500.

```
NA24631 m2           NA24631-PosCntl m2     0.0112643       0       652/1000
NA24631 m2           PTC_NA24385_S1 m2      0.0310744       0       352/1000
NA24631-PosCntl m2   PTC_NA24385_S1 m2      0.0320933       0       342/1000
```

Trying to see how it looks when comparing m2 vs m1:

```
NA24631             NA24631-PosCntl m2  0.0112643       0       652/1000  match
NA24631 m2          NA24631-PosCntl     0.0105688       0       668/1000  match
NA24631             PTC_NA24385_S1 m2   0.0370638       0       298/1000  unmatch
NA24631 m2          PTC_NA24385_S1      0.0150813       0       573/1000  unmatch
NA24631-PosCntl     PTC_NA24385_S1 m2   0.0353885       0       312/1000  unmatch
NA24631-PosCntl m2  PTC_NA24385_S1      0.0154006       0       567/1000  unmatch
```

We are seing the low concordance for mismatched samples only where PTC_NA24385 is run with m2. But this sample is low coverage in the first place - twice as lower as the other 2. So the low number might just come from that. Trying to add a higher-coverage PTC_NA24385.

```
PTC_NA24385_S11     NA24631                0.0123533       0       628/1000
PTC_NA24385_S11     NA24631-PosCntl        0.0113529       0       650/1000
PTC_NA24385_S11     PTC_NA24385_S1         0.017359        0       532/1000
PTC_NA24385_S11 m2  NA24631 m2             0.0122604       0       630/1000
PTC_NA24385_S11 m2  NA24631-PosCntl m2     0.0123533       0       628/1000
PTC_NA24385_S11 m2  PTC_NA24385_S1 m2      0.0325113       0       338/1000 
```

The high similarity between unrelated samples confirms that previously the low number came from a low PTC_NA24385 input.

### -g 3234830K

Also trying with `-g 3234830K` for pre-defined genome size to see if it affects anything:

```
time mash sketch fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz -r -g 3234830K -m 2 -o mash/NA24631-PosCntl.m2.g.msh &
Estimated genome size: 1.81814e+09
Estimated coverage:    4.27
941.35s user 1.99s system 99% cpu 15:43.34 total

time mash sketch fastqs/NA24631_S9_R1_001.fastq.gz         -r -g 3234830K -m 2 -o mash/NA24631.m2.g.msh &
Estimated genome size: 1.83433e+09
Estimated coverage:    4.025
908.76s user 2.14s system 99% cpu 15:10.91 total

time mash sketch fastqs/PTC_NA24385_S1_R1_001.fastq.gz     -r -g 3234830K -m 2 -o mash/PTC_NA24385.m2.g.msh &
```

Apparently -g option is ignored - same runtime, and same output:

```
fastqs/NA24631_S9_R1_001.fastq.gz       	fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz   0.0112643       0       652/1000
fastqs/NA24631_S9_R1_001.fastq.gz       	fastqs/PTC_NA24385_S1_R1_001.fastq.gz   	0.0310744       0       352/1000
fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz   fastqs/PTC_NA24385_S1_R1_001.fastq.gz   	0.0320933       0       342/1000
```

### Target the coverage with -c

Mash has an parameter `-c` that targets a specifiy coverage. With `-m 2 -c 1` it doesn't work, and it makes sense sinse it removes single-copy k-mers, and for a diploid genome with coverage 1x it would just zero k-mers to use. Same with `-m 2 -c 2` and `-m 2 -c 3` and `(-m 1) -c 2`. However it starts to work well with `-m 2 -c 4` or `(-m 1) -c 2`.

```
time mash sketch fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz -r -g 3234830K -c 1 -m 2 -o mash/NA24631-PosCntl.m2.g.c1.msh
Estimated genome size: 130.037
Estimated coverage:    2
Reads used:            6

time mash sketch fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz -r -g 3234830K -c 3 -m 2 -o mash/NA24631-PosCntl.m2.g.c3.msh
Estimated genome size: 1.96731e+06
Estimated coverage:    3
Reads used:            101552

time mash sketch fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz -r -g 3234830K -c 4 -m 2 -o mash/NA24631-PosCntl.m2.g.c4.msh
Estimated genome size: 1.69641e+09
Estimated coverage:    4
Reads used:            68131765
791.42s user 2.98s system 99% cpu 13:19.94 total

time mash sketch fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz -r -g 3234830K -c 2 -o mash/NA24631-PosCntl.g.c2.msh
Estimated genome size: 1.56676e+09
Estimated coverage:    2
Reads used:            28592946
347.30s user 0.71s system 99% cpu 5:48.05 total
```

With -c 2, we acheive the fastest runtime. However interesting to compare the results, so running the remaining samples:

```
time mash sketch fastqs/NA24631_S9_R1_001.fastq.gz -r -g 3234830K -c 2 -o mash/NA24631.g.c2.msh &
Estimated genome size: 1.94761e+09
Estimated coverage:    2
Reads used:            35556095
441.77s user 1.13s system 99% cpu 7:22.92 total

time mash sketch fastqs/PTC_NA24385_S1_R1_001.fastq.gz -r -g 3234830K -c 2 -o mash/PTC_NA24385.g.c2.msh &
Estimated genome size: 1.80646e+09
Estimated coverage:    1.864
Reads used:            31983302
394.25s user 0.92s system 99% cpu 6:35.21 total
```

The results are quite bad: the numbers are very close to 500/1000:

```
fastqs/NA24631_S9_R1_001.fastq.gz       	fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz   0.0202137       0       486/1000
fastqs/NA24631_S9_R1_001.fastq.gz       	fastqs/PTC_NA24385_S1_R1_001.fastq.gz   	0.0183132       0       516/1000
fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz   fastqs/PTC_NA24385_S1_R1_001.fastq.gz   	0.0214999       0       467/1000
```

-c 4 gives a runtime very close to the original run without -c, which is expected given that the estimated coverage is 4x, which is ~8x for a diploid sample. So there is no reason to experiment with that too. 

Overall, it's possible to distinguish matching samples with Mash, however it requires full 7x coverage samples to clearly do that, and won't work with a more shallow data. Also, the runtime is slower than NGSCheckMate. However, we want to see how it works with contaminati option. Sticking with NGSCheckMate for the task, but keeping mash in mind as an interesting way to estimate coverage without alignment.


### On contaminated data

```
time mash sketch fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz -r -g 3234830K -c 2 -o mash/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz.mash
# Estimated genome size: 1.94761e+09
# Estimated coverage:    2
# Reads used:            35556095
# 351.71s user 0.69s system 99% cpu 5:52.41 total

time mash sketch fastqs/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz -r -g 3234830K -c 2 -o mash/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz.mash
# Estimated genome size: 1.94761e+09
# Estimated coverage:    2
# Reads used:            35556095
# Writing to mash/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz.mash.msh...
# 354.19s user 0.51s system 99% cpu 5:54.70 total
```

Calculating the distance with `mash dist`:

```
# clean:
fastqs/NA24631_S9			fastqs/NA24631-PosCntl   0.00932572      0       698/1000
fastqs/NA24631_S9			fastqs/PTC_NA24385       0.0161637       0       553/1000
fastqs/NA24631-PosCntl   	fastqs/PTC_NA24385       0.014871        0       577/1000

# vs. 3x+0.34x contaminated
fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz (1x)       fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz   0.0202137       0       486/1000
fastqs/NA24631_S9_R1_001.fastq.gz                    fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz   0               0       1000/1000
fastqs/PTC_NA24385_S1_R1_001.fastq.gz                fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz   0.0183132       0       516/1000

# vs. 6x+0.7x contaminated

# contaminated vs contaminated
fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz   fastqs/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz    0       0       1000/1000

NA24631_S9__3x__PTC_NA24385_S11__034x  NA24631_S9
NA24631_S9__3x__PTC_NA24385_S11__034x  NA24631-PosCntl
NA24631_S9__3x__PTC_NA24385_S11__034x  PTC_NA24385_S1 (1x)
NA24631_S9__3x__PTC_NA24385_S11__034x  PTC_NA24385_S11 (3x)
NA24631_S9__3x__PTC_NA24385_S11__034x  NA24631-1KC (full)
```

Increasing k-mer size to 28:

```
time mash sketch fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz -r -p 10 -c 2 -k 28 -o mash/NA24631_S9__3x__PTC_NA24385_S11__034x_c2k28 &
time mash sketch fastqs/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz  -r -p 10 -c 2 -k 28 -o mash/NA24631_S9__6x__PTC_NA24385_S11__07x_c2k28 &
time mash sketch fastqs/NA24631_S9_R1_001.fastq.gz                     -r -p 10 -c 2 -k 28 -o mash/NA24631_S9_R1_c2k28 &  
time mash sketch fastqs/PTC_NA24385_S11_R1_001.fastq.gz                -r -p 10 -c 2 -k 28 -o mash/PTC_NA24385_S11__c2k28 &

time mash sketch fastqs/NA24631-1KC_R1.fastq.gz                        -r -p 10 -c 2 -k 28 -o mash/NA24631-1KC_R1__c2k28 &
time mash sketch fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz             -r -p 10 -c 2 -k 28 -o mash/NA24631-PosCntl_S6_R1__c2k28 &
```

That also doesn't help:

```
fastqs/NA24631-1KC_R1.fastq.gz  fastqs/NA24631-PosCntl_S6_R1_001.fastq.gz                                                       0.021441   0  378/1000
fastqs/NA24631-1KC_R1.fastq.gz  fastqs/NA24631_S9__3x__PTC_NA24385_S11__034x.fastq.gz                                           0.0188705  0  418/1000
fastqs/NA24631-1KC_R1.fastq.gz  fastqs/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz                                            0.0188705  0  418/1000
fastqs/NA24631-1KC_R1.fastq.gz  fastqs/NA24631_S9_R1_001.fastq.gz                                                               0.0188705  0  418/1000
fastqs/NA24631-1KC_R1.fastq.gz  /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/fastqs/PTC_NA24385_S11_R1_001.fastq.gz  0.020243   0  396/1000
```

Contaminated samples have same identity as clean samples. That's quite expected as they have same most abundant k-mers. So for contamination check, we can only rely on SNPs frequencies.




## ContaminationDetection tool

https://genome.sph.umich.edu/wiki/ContaminationDetection



## Conpair

Reported to perform better than VerifyBAMID and Contest.

Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048070/pdf/btw389.pdf

Github: https://github.com/nygenome/conpair

Dependencies: "GATK 2.3 or higher.". Breaks with GATK 4. Here we go again...

```
cd /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/conpair
source load.sh

./conpair/scripts/run_gatk_pileup_for_sample.py -B /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bams/NA24631_S9__3x__PTC_NA24385_S11__034x-sort.bam -O NA24631_S9__3x__PTC_NA24385_S11__034x.pileup

./conpair/scripts/run_gatk_pileup_for_sample.py -B /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/bams/NA24631-1KC-ready.bam -O NA24631-1KC.pileup

python ./conpair/scripts/verify_concordance.py -T NA24631-1KC.pileup -N NA24631_S9__3x__PTC_NA24385_S11__034x.pileup
0.844
Based on 1732/7387 markers (coverage per marker threshold: 10 reads)
Minimum mappinq quality: 10
Minimum base quality: 20
```

Gives 0.844. Seems legit.





