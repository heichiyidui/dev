#!/usr/bin/env python3
import sys

ifile=open('index/cath_s35.seq')
seqs={}
for line in ifile:
    dom_id = line[1:-1]
    line = ifile.readline()[:-1]
    seqs[dom_id] = line
ifile.close()

for dom_id in seqs.keys():
    dssp_ss  = ['X'] * len(seqs[dom_id])
    dssp_acc = ['NA'] * len(seqs[dom_id])

    ifile=open('dssp_ss/'+dom_id)
    for line in ifile:
        if line.startswith('  #  RESIDUE AA'):
            break;

    for line in ifile:
        if line[6:10] =='    ':
            continue

        res_num = int(line[6:10]) - 1

        res_AA = line[13]
        if res_AA != seqs[dom_id][res_num]:
            print('What? Sequence mismatching at:',dom_id,line[:-1])

        ss = line[16]
        if ss == ' ':   ss='C'
        dssp_ss[res_num]=ss

        acc = line[35:38]
        dssp_acc[res_num] = acc

    ifile.close()

    print('>'+dom_id )
    # print(''.join(dssp_ss))
    # print(' '.join(dssp_acc))