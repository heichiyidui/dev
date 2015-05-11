#!/usr/bin/env python3 
import sys

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

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

    #########################
    # read DSSP file 

    dssp_acc=['NA']* len(cseq)

    ifile=open('dssp_ss/'+dom_id)
    for line in ifile:
        if line.startswith('  #  RESIDUE AA'):
            break;
    
    for line in ifile:
        if line[6:10] =='    ':
            continue
        res_num=int(line[6:10])-1
        acc=line[35:38]
        dssp_acc[dom_num_2_seq_num[res_num]]=acc 
    ifile.close()

    print('>'+dom_id)
    print(' '.join(dssp_acc))
