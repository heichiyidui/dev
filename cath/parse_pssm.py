#!/usr/bin/env python3 

domain_ids=open('t.ls').read().split()

for dom_id in domain_ids:
    #########################
    # read consensus sequence
    ifile = open('cseq/'+dom_id)
    cseq=''.join(ifile.read().split('\n')[1:])
    ifile.close()

    Cseq=cseq.upper()

    #########################
    # read PSSM file 
    pssm=[]

    ifile=open('pssm/'+dom_id)
    for line in ifile:
        if line.startswith('           A  R  N  D  C  Q  E  G'):
            break;
    
    for line in ifile:
        if line=='\n':
            break;
        res_num=int(line[0:5])-1
        aa = line[6]
        if Cseq[res_num] != aa:
            print(dom_id,res_num)
        pssm_aa=[]
        for i in range(9,69,3):
            pssm_aa.append(line[i:i+3].strip())
        pssm.append('|'.join(pssm_aa))
    ifile.close()

    print(">"+dom_id)
    print(' '.join(pssm))

