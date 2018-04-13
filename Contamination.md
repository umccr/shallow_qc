Contamination control from shallow FASTQ
========================================

Out idea is to understand if we can perform a contamination pre-flight check from shallow data at 3x-7x coverage. We will have shallow GiaB samples data with an unknown contamination, and we can compare it to a full GiaB data known in advance. There are 2 basic approaches we want to explore: 

1. Call variants and compare AFs against a truth. Contaminated data should result in a consistently lower AF of each score, e.g. for 10% contamination, AFs should be on avarage 10% lower.
2. Count k-mers in FASTQ, the k-mer identity should be 10% less than expected for a clean sample.

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
time gunzip -c NA24631_S9_R1_001.fastq.gz      | head -n 66699966 | gzip -c > NA24631_S9_R1_001.head_3x.fastq.gz                # 2m35s
time gunzip -c PTC_NA24385_S11_R1_001.fastq.gz | head -n 7411107  | gzip -c > PTC_NA24385_S11_R1_001.head_034x.fastq.gz
```

The coverage might be quite low, so we will also test 6x+0.7x (129393200 + 15095873) by using the R2 reads as well:

```
# 6x + 0.7x:
time gunzip -c NA24631_S9_R1_001.fastq.gz      | head -n 66699966 | gzip -c >  NA24631_S9.head_6x.fastq.gz
time gunzip -c NA24631_S9_R2_001.fastq.gz      | head -n 66699966 | gzip -c >> NA24631_S9.head_6x.fastq.gz
time gunzip -c PTC_NA24385_S11_R1_001.fastq.gz | head -n 7411107  | gzip -c >  PTC_NA24385_S11.head_07x.fastq.gz
time gunzip -c PTC_NA24385_S11_R2_001.fastq.gz | head -n 7411107  | gzip -c >> PTC_NA24385_S11.head_07x.fastq.gz
```

For the sake of experiment, aslo trying to subsample correctly just in case:

```
time seqtk sample -s 11 PTC_NA24385_S11_R1_001.fastq.gz 0.9 | gzip -c > PTC_NA24385_S11_R1_001.seqtk_90pct.fastq.gz
# seqtk sample -s 11 PTC_NA24385_S11_R1_001.fastq.gz 0.9  	343.89s user 20.75s system 10% cpu 1:00:30.13 total
# gzip -c > PTC_NA24385_S11_R1_001.seqtk_90pct.fastq.gz  	3100.37s user 13.22s system 85% cpu 1:00:30.15 total
# Timing for seqtk only: 
# seqtk sample -s 11 PTC_NA24385_S11_R1_001.fastq.gz 0.9    330.30s user 17.57s system 61% cpu 9:24.57 total
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
# 920.55s user 10.03s system 97% cpu 15:49.92 total

time /data/cephfs/punim0010/projects/Saveliev_Fingerprinting/NGSCheckMate/ngscheckmate_fastq -p 9 -1 fastqs/NA24631_S9__6x__PTC_NA24385_S11__07x.fastq.gz NGSCheckMate/SNP/SNP.pt > ncm_files/NA24631_S9__6x__PTC_NA24385_S11__07x.ncm
Number of reads in the file : 37055537
1514.89s user 25.62s system 97% cpu 26:23.89 total

python2 NGSCheckMate/vaf_ncm.py -f -I ncm_files -O output_contaminated

# vs clean B
NA24631-PosCntl_R1.fastq.gz            matched    NA24631_S9_R2_001.fastq.gz             0.8014  1.42   # A - 3.8x  vs.  B - 3.6x
NA24631-PosCntl_D1_R2_001.fastq.gz     matched    NA24631_S9_R2_001.fastq.gz             0.7441  0.43   # A - 1x    vs.  B - 3.6x 
PTC_NA24385_S1_R2_001.fastq.gz         unmatched  NA24631_S9_R2_001.fastq.gz             0.1174  0.67   # C - 1.6x  vs.  B - 3.6x
PTC_NA24385_S1_SS025_R2_001.fastq.gz   unmatched  NA24631_S9_R2_001.fastq.gz             0.103   0.2    # C - 1x    vs.  B - 3.6x

# vs 3x+0.34x contaminated
NA24631-PosCntl_R1.fastq.gz            matched    NA24631_S9__3x__PTC_NA24385_S11__034x  0.7272  0.37   # A - 3.8x  vs.  B - 3x + C - 0.34x
NA24631-PosCntl_D1_R2_001.fastq.gz     matched    NA24631_S9__3x__PTC_NA24385_S11__034x  0.6691  0.37   # A - 1x    vs.  B - 3x + C - 0.34x
NA24631_S9_R2_001.fastq.gz             matched    NA24631_S9__3x__PTC_NA24385_S11__034x  0.8778  0.37   # B - 3.6x  vs.  B - 3x + C - 0.34x
NA24631_S9_SS025_R2_001.fastq.gz       matched    NA24631_S9__3x__PTC_NA24385_S11__034x  0.7794  0.37   # B - 1x    vs.  B - 3x + C - 0.34x
PTC_NA24385_S1_R2_001.fastq.gz         unmatched  NA24631_S9__3x__PTC_NA24385_S11__034x  0.0974  0.37   # C - 1.6x  vs.  B - 3x + C - 0.34x  
PTC_NA24385_S1_SS025_R2_001.fastq.gz   unmatched  NA24631_S9__3x__PTC_NA24385_S11__034x  0.0961  0.2    # C - 1x    vs.  B - 3x + C - 0.34x

# vs 6x+0.7x contaminated
NA24631-PosCntl_R1.fastq.gz            matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.6793  0.4    # A - 3.8x  vs.  B - 6x + C - 0.7x
NA24631-PosCntl_D1_R2_001.fastq.gz     matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.6224  0.4    # A - 1x    vs.  B - 6x + C - 0.7x
NA24631_S9_R2_001.fastq.gz             matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.816   0.4    # B - 3.6x  vs.  B - 6x + C - 0.7x
NA24631_S9_SS025_R2_001.fastq.gz       matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.7357  0.4    # B - 1x    vs.  B - 6x + C - 0.7x
PTC_NA24385_S1_R2_001.fastq.gz         unmatched  NA24631_S9__6x__PTC_NA24385_S11__07x   0.1488  0.4    # C - 1.6x  vs.  B - 6x + C - 0.7x
PTC_NA24385_S1_SS025_R2_001.fastq.gz   unmatched  NA24631_S9__6x__PTC_NA24385_S11__07x   0.1403  0.2    # C - 1x    vs.  B - 6x + C - 0.7x

# contaminated vs contaminated
NA24631_S9__3x__PTC_NA24385_S11__034x  matched    NA24631_S9__6x__PTC_NA24385_S11__07x   0.9836  0.37   # B - 3x + C - 0.34x  vs B - 6x + C - 0.7x
```

We see that 6x+0.7x yields a higher similarity with PTC_NA24385 (~0.14) than when comparing against a pure NA24631 (~0.11), and at the same time it yields a lower similarity with NA24631 (~0.6) than when comparing against a pure NA24631 (~0.8) - which looks promising, however doesn't give us exact thresholds to confidently tell one from the other. We will dig into the VAF files to see if we can extract the contamination rate from them - e.g., calculate a median differece between AFs for all matching calls.




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
fastqs/NA24631_S9			fastqs/NA24631-PosCntl	 0.00932572      0       698/1000
fastqs/NA24631_S9			fastqs/PTC_NA24385	     0.0161637       0       553/1000
fastqs/NA24631-PosCntl_S6	fastqs/PTC_NA24385	     0.014871        0       577/1000
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
NA24631 m2       	 NA24631-PosCntl m2     0.0112643       0       652/1000
NA24631 m2     	     PTC_NA24385_S1 m2 	    0.0310744       0       352/1000
NA24631-PosCntl m2   PTC_NA24385_S1 m2 		0.0320933       0       342/1000
```

Trying to see how it looks when comparing m2 vs m1:

```
NA24631       		NA24631-PosCntl m2  0.0112643       0       652/1000  match
NA24631 m2     		NA24631-PosCntl     0.0105688       0       668/1000  match
NA24631       		PTC_NA24385_S1 m2   0.0370638       0       298/1000  unmatch
NA24631 m2			PTC_NA24385_S1		0.0150813       0       573/1000  unmatch
NA24631-PosCntl	    PTC_NA24385_S1 m2   0.0353885       0       312/1000  unmatch
NA24631-PosCntl m2  PTC_NA24385_S1		0.0154006       0       567/1000  unmatch
```

We are seing the low concordance for mismatched samples only where PTC_NA24385 is run with m2. But this sample is low coverage in the first place - twice as lower as the other 2. So the low number might just come from that. Trying to add a higher-coverage PTC_NA24385.

```
time mash sketch fastqs/PTC_NA24385_S11_R1_001.fastq.gz -r -p 4 -m 2 -o mash/PTC_NA24385_S11_R1_001.m2.fastq.gz & 
time mash sketch fastqs/PTC_NA24385_S11_R1_001.fastq.gz -r -p 4 -o mash/PTC_NA24385_S11_R1_001.fastq.gz & 

PTC_NA24385_S11		NA24631			       0.0123533       0       628/1000
PTC_NA24385_S11		NA24631-PosCntl	       0.0113529       0       650/1000
PTC_NA24385_S11		PTC_NA24385_S1	       0.017359        0       532/1000
PTC_NA24385_S11	m2	NA24631 m2 			   0.0122604       0       630/1000
PTC_NA24385_S11	m2	NA24631-PosCntl m2	   0.0123533       0       628/1000
PTC_NA24385_S11	m2	PTC_NA24385_S1 m2 	   0.0325113       0       338/1000
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

Overall, it's possible to distinguish matching samples with Mash, however it requires full 7x coverage samples to clearly do that, and won't work with a more shallow data, so we can't it expect to work with contaminated data. Also, the runtime is slower, so NGSCheckMate looks like a better option. Sticking with NGSCheckMate for the task, but keeping mash in mind as an interesting way to estimate coverage without alignment.







