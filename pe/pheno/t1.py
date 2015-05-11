#! /usr/bin/env python3

ifile=open('data/Quant_traits_resid.csv')
pheno={}
ifile.readline()
for line in ifile:
    if line.split(',')[30] != 'NA':
        pheno[line.split(',')[0]]=line.split(',')[30]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()
