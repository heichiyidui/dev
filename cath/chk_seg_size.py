#!/usr/bin/env python3 

import sys
import numpy

domain_ids=open(sys.argv[1]).read().split()

for id in domain_ids:
    ifile=open('dompdb/'+id)
    atom_lines=ifile.read().split('\n')[:-1]
    ifile.close()
    
    res_C=[]
    res_N=[]
    for line in atom_lines:
        if line[12:16]==' C  ':
            res_C.append( numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
        if line[12:16]==' N  ':
            res_N.append( numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
    seg_lengths=[]
    seg_len = 1
    for i in range(1,len(res_C)):
        if numpy.linalg.norm(res_C[i-1]-res_N[i]) > 2.5:
            seg_lengths.append(seg_len)
            seg_len = 1
        seg_len += 1
    seg_lengths.append(seg_len)
    
    print(id, len(res_C), len(seg_lengths), end= ' ')
    for seg_len in seg_lengths:
        print(seg_len, end =' ')
    print()
