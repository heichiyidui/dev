#!/usr/bin/env python3 
import sys
import numpy

dom_ids=open('t.ls').read().split()

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
    res_CA=[]
    
    for line in atom_lines:
        if line[12:16] == ' C  ':
            res_C.append( numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
        if line[12:16] == ' N  ':
            res_N.append( numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
        
        if line[12:16] == ' CA ':
            res_CA.append(numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
        
    
    breaks=set([])
    for i in range(1,len(res_ids)):
        if numpy.linalg.norm(res_C[i-1]-res_N[i]) > 2.5:
            breaks.add(i-1)           
    
    for i in range(len(res_ids)-3):
        if i in breaks:
            continue
        if i+1 in breaks:
            continue
        if i+2 in breaks:
            continue
        print( numpy.linalg.norm(res_CA[i  ]-res_CA[i+1]),\
               numpy.linalg.norm(res_CA[i+1]-res_CA[i+2]),\
               numpy.linalg.norm(res_CA[i+2]-res_CA[i+3]))
