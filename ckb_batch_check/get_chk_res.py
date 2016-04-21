#!/usr/bin/env python3

good_snps = open('t.in').read().split()
ifile=open('to_exam.ls')
for line in ifile:
    snp_id = line[:-1]
    if snp_id not in good_snps:
        print(snp_id,0)
    else:
        print(snp_id,1)
