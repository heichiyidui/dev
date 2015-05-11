#!/usr/bin/env python3 
import sys
import numpy

dom_ids=[]
ifile=open(sys.argv[1])
for line in ifile:
    cols=line.split()
    dom_ids.append(cols[0])
ifile.close()

for dom_id in dom_ids:
    ifile=open('../cath/dompdb/'+dom_id)
    atom_lines=ifile.read().split('\n')[:-1]
    ifile.close()
    
    res_ids=[]
    for line in atom_lines:
        if line[21:27] not in res_ids:
            res_ids.append(line[21:27])
    
    res_C=[]
    res_N=[]
    for line in atom_lines:
        if line[12:16] == ' C  ':
            res_C.append( numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
        if line[12:16] == ' N  ':
            res_N.append( numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
    
    seg_lens=[]
    seg_len=1
    for i in range(1,len(res_ids)):
        if numpy.linalg.norm(res_C[i-1]-res_N[i]) > 2.5:
            seg_lens.append(seg_len)            
            seg_len=1
            continue
        seg_len += 1
    seg_lens.append(seg_len)
    
    print(dom_id, '{:4d}'.format(len(res_ids)), len(seg_lens), end=' ')
    
    for seg_len in seg_lens:
        print('{:4d}'.format(seg_len),end = ' ')
    print()
