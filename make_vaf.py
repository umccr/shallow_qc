import sys

txt = '''9       90500405        .       A       C,<*>   0       .       DP=62;I16=0,0,26,36,0,0,2280,84038,0,0,3720,223200,0,0,1366,31882;QS=0,1,0;VDB=0.129645;SGB=-0.693147;MQSB=1;MQ0F=0     PL:DP:ADF:AD    255,187,0,255,187,255:62:0,26,0:0,62,0
NA NA24631_S9__3x__PTC_NA24385_S11__034x-sort.mpileup_10:100018844-100018844'''

'''
# making this:
index   ref     alt     vaf
0       0       0       NA
1       NA      NA      NA
2       0       2       1.000000
'''

if len(sys.argv) > 1:
    with open(sys.argv[1]) as f:
        txt = f.read()

record_by_coord = dict()

print('index\tref\talt\tvaf')

for line in txt.split('\n'):
    if line.strip():
        fs = line.split()
        if fs[0] == 'NA':
            coord = fs[1].split('_')[-1]
            chrom, pos = coord.split(':')
            pos = int(pos.split('-')[0])
            record_by_coord[(chrom, pos)] = ['NA', 'NA', 'NA']
        else:
            chrom, pos = fs[0], int(fs[1])
            ad = fs[9].split(':')[-1]
            ad_ref, ad_alt = ad.split(',')[0:2]
            summ = int(ad_ref) + int(ad_alt)
            record_by_coord[(chrom, pos)] = [ad_ref, ad_alt, (int(ad_alt) / summ) if summ > 0 else 0]


for i, (coord, rec) in enumerate(sorted(record_by_coord.items())):
    print('\t'.join(map(str, [i] + rec)))

