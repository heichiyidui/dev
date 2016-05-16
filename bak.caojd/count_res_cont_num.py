#!/usr/bin/env python3
import sys
while '/share/python/lib64/python' in sys.path :
    sys.path.remove('/share/python/lib64/python') 
# the cluster got this wrong dir

import numpy as np 
res_cont_sums=[0]*15

dom_ids=open('index/s.ls').read().split()

ifile=open('index/cath_s35.condef')
for line in ifile:
    c_line=ifile.readline().strip()
    
    dom_id = line.split()[0][1:]
    if dom_id not in dom_ids: continue 
    
    dom_len = int (line.split()[1])
    cont_num = int (line.split()[2])
    cols=c_line.split()
    
    cont_map=np.zeros((dom_len,dom_len))
    
    for i in range(cont_num):
        c_i=int(cols[i*2])
        c_j=int(cols[i*2+1])
        cont_map[c_i,c_j]=1
        cont_map[c_j,c_i]=1
    
    for s in sum(cont_map):
        res_cont_sums[int(s)]+=1

ifile.close()    

for i in range(len(res_cont_sums)):
    print(i,res_cont_sums[i])

print()

print(sum(res_cont_sums))

print()
t_sum=0
for i in range(len(res_cont_sums)):
    t_sum += res_cont_sums[i]
    print(i,t_sum)
