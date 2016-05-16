#!/usr/bin/env python3
import collections

AA_CODE="ARNDCQEGHILKMFPSTWYV"
AA_to_INT = {}
for i in range(20):
    AA_to_INT[AA_CODE[i]] = i

AA_AA_to_INT = {}
for i in range(20):
    for j in range(20):
        AA_AA_to_INT[ AA_CODE[i],AA_CODE[j] ] = i*20 +j

#######################################
# read the sequences
dom_seq={}
ifile=open('../cath/index/cath_s35.seq')
for line in ifile:
    dom_id = line[1:-1]
    seq = ifile.readline()[:-1]
    dom_seq[dom_id] = seq
ifile.close()

#######################################
# read the contact matrices
cont_mats = {}
ifile=open('../cath/index/cath_s35.condef')
for line in ifile:
    dom_id = line[1:-1]
    cont_mat = set([])
    cols = ifile.readline()[:-1].split()
    for col in cols:
        c_i,c_j = col.split('-')
        i,j = int(c_i),int(c_j)
        cont_mat.add((i,j))
    cont_mats[dom_id] = cont_mat
ifile.close()

#######################################
# count the contacts

aa_aa_counter = collections.Counter()

dom_ls = open('index/dom.ls').read().split()

for dom_id in dom_ls:
    conts=[]
    seq = dom_seq[dom_id]
    for (i,j) in cont_mats[dom_id]:
        if j-i < 5 : continue
        conts.append(AA_AA_to_INT[seq[i],seq[j]])
    aa_aa_counter.update(conts)

sum_aa_aa = sum(aa_aa_counter.values())

for c1 in AA_CODE:
    for c2 in AA_CODE:
        freq = aa_aa_counter[AA_AA_to_INT[c1,c2]]/sum_aa_aa
        print('{:0.8f}'.format(freq),end=' ')
    print()


