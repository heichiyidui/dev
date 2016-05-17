#!/usr/bin/env python3

MAX_DIS =  253

#######################################
# maps
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
        if j-i < 5:
            continue
        cont_mat.add((i,j))
    cont_mats[dom_id] = cont_mat
ifile.close()

#######################################
# read alignment distances
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
# sum the alignments
import numpy as np
Pi = np.genfromtxt('cao/Pi_5')
Pi.shape = (400)

sum_aln_mats=[]
for i in range(100):
    sum_aln_mats.append(np.outer(Pi,Pi))    # Laplace smoothing?

dom_ls = open('index/dom.ls').read().split()

for dom_id in dom_ls:
    print('>'+dom_id)
    ifile=open('../cath/bl_out/'+dom_id)
    ifile.readline()
    dom_seq = ifile.readline()[:-1]
    for dis in dom_dis[dom_id]:
        line = ifile.readline()
        aln_seq = ifile.readline()[:-1]
        dis_index = int(dis/(MAX_DIS/100))
        print(dis_index)
        for (i,j) in cont_mats[dom_id]:
            if aln_seq[i] == '-':
                continue
            if aln_seq[j] == '-':
                continue
            dom_cont = AA_AA_to_INT[dom_seq[i],dom_seq[j]]
            aln_cont = AA_AA_to_INT[aln_seq[i],aln_seq[j]]
            sum_aln_mats[dis_index][dom_cont][aln_cont] += 1
            sum_aln_mats[dis_index][aln_cont][dom_cont] += 1
    ifile.close()

#######################################
# normalize the sum matrices, write them to files
for i in range(100):
    sum_aln_mats[i] /= np.sum(sum_aln_mats[i])
    out_file_name = 'sum_mat/{:03d}'.format(i)
    np.savetxt(out_file_name,sum_aln_mats[i],'%08e')
