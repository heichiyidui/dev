#!/usr/bin/env python3
import sys

AA_NAME_TO_CODE={'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',\
                 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',\
                 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',\
                 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

ifile=open('index/cath_s35.seq')
seqs={}
for line in ifile:
    dom_id = line[1:-1]
    line = ifile.readline()[:-1]
    seqs[dom_id] = line
ifile.close()

for dom_id in seqs.keys():
    stride_ss  = ['X']  * len(seqs[dom_id])
    stride_acc = ['NA'] * len(seqs[dom_id])

    ifile=open('stride_ss/'+dom_id)
    for line in ifile:
        if not line.startswith("ASG"):
            continue
        code=AA_NAME_TO_CODE[line[5:8]]
        res_num=int(line[16:20])-1

        if code != seqs[dom_id][res_num]:
            print('What? Sequence mismatching at:',dom_id,line[:-1])

        ss = line[24]
        if ss == 'b':  ss = 'B'

        stride_ss[res_num] = ss

        acc = line[64:69]
        stride_acc[res_num] = acc
    ifile.close()

    print(">"+dom_id)
    # print(''.join(stride_ss))
    # print(' '.join(stride_acc))