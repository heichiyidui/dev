#!/usr/bin/env python3 


pheno={}
ifile=open('input_files/ba3.pheno')
for line in ifile:
    cols=line[:-1].split()
    pheno[cols[0]]=cols[1]
ifile.close()

AD_ids=[]
CTL_ids=[]
for s_id in pheno.keys():
    if pheno[s_id] == 'AD':
        AD_ids.append(s_id)
    if pheno[s_id] == 'CTL':
        CTL_ids.append(s_id)

import random
random.seed(514)
random.shuffle(AD_ids)

AD_ids=AD_ids[:59]
# finished reading AD and CTL ids 

import numpy

ifile=open('input_files/ba3_50_probes.txt')
line=ifile.readline()
print(line[:-1])
subject_ids=line[:-1].split()[1:]

for line in ifile:
    cols=line[:-1].split()
    probe_id=cols[0]
    
    exprs_profile=[]
    for col in cols[1:]:
        exprs_profile.append(float(col))
    
    exprs_profile_subset=[]
    for i in range(len(exprs_profile)):
        if subject_ids[i] in AD_ids or subject_ids[i] in CTL_ids:
            exprs_profile_subset.append(exprs_profile[i])
    
    exprs_mean = numpy.average(exprs_profile_subset)
    exprs_std  = numpy.std(exprs_profile_subset)
    
    print(probe_id,end=' ')
    for i in range(len(exprs_profile)):
        print((exprs_profile[i] - exprs_mean ) / exprs_std, end=' ')
    print()

ifile.close()
    