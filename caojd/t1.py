#!/usr/bin/env python3

#######################################
# 1. read alignment distances
dom_dis={}
ifile=open('dis.out')
line = ifile.readline()
dom_id = line[1:-1]
while True:
    line = ifile.readline()
    if not line : break;
    if line.startswith('>'):
        dom_id = line[1:-1]
        continue
    if dom_id not in dom_dis.keys():
        dom_dis[dom_id]=[]
    dom_dis[dom_id].append(float(line[:-1]))


#######################################
# 2. read alignments, sum up the aln matrices

import numpy as np

sum_aln_mats=[]
for i in range(100):
    sum_mat=[]
    for j in range(20):
        sum_mat.append([0]*20)
    sum_mat = np.array(sum_mat)

    sum_aln_mats.append(sum_mat)


