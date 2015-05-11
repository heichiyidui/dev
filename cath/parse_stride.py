#!/usr/bin/env python3 
import sys

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

AA_NAME_TO_CODE={'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',\
                 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',\
                 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',\
                 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

for dom_id in domain_ids:
    #########################
    # read consensus sequence
    ifile = open('cseq/'+dom_id)
    cseq=''.join(ifile.read().split('\n')[1:])
    ifile.close()

    dom_num_2_seq_num={}
    j=0
    for i in range(len(cseq)):
        if cseq[i].isupper():
            dom_num_2_seq_num[j]=i
            j+=1

    dseq=''
    for c in cseq:
        if c.isupper():
            dseq += c 
    #########################
    # read STRIDE file 

    stride_ss=['X']* len(cseq)

    ifile=open('stride_ss/'+dom_id)
    for line in ifile:
        if not line.startswith("ASG"):
            continue
        code=AA_NAME_TO_CODE[line[5:8]]
        res_num=int(line[16:20])-1
        # if code != dseq[res_num]:
        #     print(id,res_num)
        ss=line[24]
        if ss=='b':ss='B'
        stride_ss[dom_num_2_seq_num[res_num]]=ss 
    ifile.close()

    print(">"+dom_id)
    print(''.join(stride_ss))