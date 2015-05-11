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

    dseq=''
    for c in cseq:
        if c.isupper():
            dseq += c 
    #########################
    # read contact file 
    contacts=[]
    ifile=open('conDef/'+dom_id)
    for line in ifile:
        contc=line.strip().replace(' ','-')
        contacts.append(contc)
    ifile.close()

    print(">"+dom_id)
    print(' '.join(contacts))
