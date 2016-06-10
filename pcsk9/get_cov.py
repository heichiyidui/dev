#!/usr/bin/env python3

################################################################################
# 1. read the FID, IID, rc, age, sex from pheno.csv                            #
# 2. read the 10 PCs from plink.eigenvec                                       #
# 3. merge and output the table                                                #
################################################################################

#######################################
# 1. read the FID, IID, rc, age, sex from pheno.csv
rc_2_int={}
rc_2_int['12']=0
rc_2_int['16']=1
rc_2_int['26']=2
rc_2_int['36']=3
rc_2_int['46']=4
rc_2_int['52']=5
rc_2_int['58']=6
rc_2_int['68']=7
rc_2_int['78']=8
rc_2_int['88']=9
FID={}
IID=[]
rc={}
age={}
sex={}

ifile=open('pheno.csv')
ifile.readline()
for line in ifile:
    cols = line[:-1].split('\t')
    iid      = cols[1]
    IID.append(iid)
    FID[iid] = cols[0]
    rc [iid] = rc_2_int[cols[3]]
    age[iid] = cols[5]
    sex[iid] = cols[6]

ifile.close()

#######################################
# 2. read the 10 PCs
ifile=open('plink.eigenvec')
pcs={}
for line in ifile:
    cols = line[:-1].split()
    iid=cols[1]
    pcs[iid]=cols[2:]

ifile.close()

#######################################
# 3. merge and output
print('FID IID age sex',end = ' ')
print('rc1 rc2 rc3 rc4 rc5 rc6 rc7 rc8 rc9',end=' ')
print('pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10',end='\n')
for iid in IID:
    print(FID[iid],iid,age[iid],sex[iid],end=' ')
    rcs=[0]*9
    if rc[iid] != 9:
        rcs[rc[iid]]=1
    print(' '.join([str(x) for x in rcs]), end =' ')
    print(' '.join(pcs[iid]),end = '\n')