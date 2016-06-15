#!/usr/bin/env python3
import sys

# read SNP allele 1 and 2
ifile=open('ckb_ph12_s3.bim')
SNP_A1={}
SNP_A2={}
for line in ifile:
    cols = line[:-1].split()
    SNP_A1[cols[1]] = cols[4]
    SNP_A2[cols[1]] = cols[5]
ifile.close()

# insert the A2 column
ifile=open(sys.argv[1]) # plink assoc linear file
cols=ifile.readline()[:-1].split()

is_A2_found=False
if 'A2' in cols:
    sys.stderr.write('A2 is there already')
    is_A2_found=True

if not is_A2_found:
    cols.insert(4,'A2')
print(' '.join(cols))
for line in ifile:
    cols = line[:-1].split()
    snp_id = cols[1]
    if not is_A2_found:
        if SNP_A1[snp_id] == cols[3]:
            cols.insert(4,SNP_A2[snp_id])
        else:
            cols.insert(4,SNP_A1[snp_id])
    print(' '.join(cols))
ifile.close()