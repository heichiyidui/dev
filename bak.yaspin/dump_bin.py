#!/usr/bin/env python3
import pickle

ifile=open('index/domain.ls')
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for domain_id in domain_ids:
    ifile=open('in/'+domain_id)
    pssm=[]
    ss_def=[]
    for line in ifile:
        ss_def.append(line[4])
        for i in range(20):
            pssm.append(float(line[11+i*3:14+i*3]))
    ifile.close()
    
    ofile=open('bin/'+domain_id,'wb')
    pickle.dump(pssm,ofile)
    pickle.dump(ss_def,ofile)
    ofile.close()
    

