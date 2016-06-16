#!/usr/bin/env python3

ifile=open('plink.map')
SNPs=[]
for line in ifile:
    cols = line[:-1].split()
    SNPs.append(cols[1])
ifile.close()

print('FID IID',end=' ')
print(' '.join(SNPs))

ifile=open('plink.ped')
for line in ifile:
    cols = line[:-1].split()
    calls = [int(x) for x in cols[6:]]

    dosage=[0]*len(SNPs)
    for i in range(len(SNPs)):
        dosage[i] = calls[i*2] + calls[i*2+1] - 2
        if 0 in [calls[i*2] , calls[i*2+1] ]:
            dosage[i] = -9

    str_dosage = [str(x) for x in dosage]
    for n,i in enumerate(str_dosage):
        if i=='-9':
            str_dosage[n]='NA'

    print(cols[0],cols[1],end=' ')
    print(' '.join(str_dosage))
ifile.close()
