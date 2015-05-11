#!/usr/bin/env python3

dssp_defs={}
ifile=open('index/cath_s35.dssp.acc')
for line in ifile:
    dom_id=line.strip()[1:]
    line=ifile.readline().strip()
    dssp_defs[dom_id]=line 
ifile.close()

stride_defs={}
ifile=open('index/cath_s35.stride.acc')
for line in ifile:
    dom_id=line.strip()[1:]
    line=ifile.readline().strip()
    stride_defs[dom_id]=line 
ifile.close()

import numpy
for dom_id in dssp_defs.keys():
    dssp_acc=[]
    for acc_str in dssp_defs[dom_id].split():
        dssp_acc.append(float(acc_str))
    
    stride_acc=[]
    for acc_str in stride_defs[dom_id].split():
        stride_acc.append(float(acc_str))
    
    print(dom_id,numpy.corrcoef(dssp_acc,stride_acc)[1,0])
