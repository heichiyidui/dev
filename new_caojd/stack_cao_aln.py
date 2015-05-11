#!/usr/bin/env python3 
##!/home/klinbrc/bin/python3 
#

MAX_ALN_DEPTH  = 2000

AA_CODE_TO_INT = {'A': 0 ,'R': 1 ,'N': 2 ,'D': 3 ,'C': 4 ,\
                  'Q': 5 ,'E': 6 ,'G': 7 ,'H': 8 ,'I': 9 ,\
                  'L': 10,'K': 11,'M': 12,'F': 13,'P': 14,\
                  'S': 15,'T': 16,'W': 17,'Y': 18,'V': 19,};

import sys

dom_ids=open(sys.argv[1]).read().split()

# read alignment depths and domain contact definitions
dom_depths={}
ifile=open('index/cath_s35.aln_depth')
for line in ifile:
    cols=line.split()
    if cols[0] not in dom_ids:
        continue
    dom_depths[cols[0]]=int(cols[1])
ifile.close()

dom_contacts={}
ifile=open('index/cath_s35.condef')
for line in ifile:
    dom_id=line.strip()[1:]
    cont_line=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    cols=cont_line.split()
    contacts=[]
    for col in cols:
        (c_i,c_j)=(int(col.split('-')[0]),int(col.split('-')[1]))
        contacts.append((c_i,c_j))
    dom_contacts[dom_id]=contacts
ifile.close()

# read distance thresholds
ifile=open('t_dis_split.ls')
dis_split_list=[]
for line in ifile:
    dis_split_list.append(float(line.strip()))
ifile.close() # split alignments into 100 different sum matrices 

sum_aln_mats=[]
for i in range(100):
    sum_aln_mats.append([0]*160000)

import os
import subprocess
import random
random.seed(514)

# for each domain 
for dom_id in dom_ids:
    
    # 1 read alignment 
    # 2 put vtml mat to a file 
    # 3 get vtml dis
    # 4 add cao mat to the sum 
    
    ifile=open('compact_paln/'+dom_id)
    ifile.readline()
    dom_seq=ifile.readline().strip()
    dom_len=len(dom_seq)
    
    for line in ifile:
        subjc_seq=ifile.readline().strip()
        
        if dom_depths[dom_id] > MAX_ALN_DEPTH:
            if random.random() > MAX_ALN_DEPTH/dom_depths[dom_id]:
                continue
        
        vtml_aln=[0]*400
        for i in range(dom_len):
            if subjc_seq[i] == '-':
                continue
            dom_i = AA_CODE_TO_INT[dom_seq[i]]
            sub_i = AA_CODE_TO_INT[subjc_seq[i]]
            
            vtml_aln[dom_i*20+sub_i] += 1
            vtml_aln[sub_i*20+dom_i] += 1
        
        ofile=open(sys.argv[1]+'.vmat','w')
        ofile.write(' '.join(str(x) for x in vtml_aln)+'\n')
        ofile.close()
        
        dis_out=subprocess.check_output('/usr/bin/Rscript vmat_dis.R ' \
                                        +sys.argv[1]+'.vmat',shell=True)
        # dis_out=subprocess.check_output('/share/apps/R_3.0.2/bin/Rscript' + \
        #                 ' vmat_dis.R '+ sys.argv[1]+'.vmat',shell=True)
        vtml_dis=float(dis_out.split()[1])
        # got VTML distance
        
        if vtml_dis < 1.5:
            continue 
        
        aln_class=0
        for i in range(len(dis_split_list)):
            if vtml_dis > dis_split_list[i]:
                aln_class += 1
        
        for contact in dom_contacts[dom_id]:
            c_i=contact[0]
            c_j=contact[1]
            # luckly global contacts are not affected by chain breaks
            if c_j - c_i < 5 : 
                continue
            if subjc_seq[c_i]=='-':
                continue
            if subjc_seq[c_j]=='-':
                continue
            cont_dom   = AA_CODE_TO_INT[dom_seq[c_i]]*20 \
                        +AA_CODE_TO_INT[dom_seq[c_j]]
            
            cont_subjc = AA_CODE_TO_INT[subjc_seq[c_i]]*20 \
                        +AA_CODE_TO_INT[subjc_seq[c_j]]
            
            sum_aln_mats[aln_class][cont_dom*400 + cont_subjc] += 1
            sum_aln_mats[aln_class][cont_subjc*400 + cont_dom] += 1
        
    ifile.close()

# output sum 
for i in range(100):
    for j in range(160000):
        print(sum_aln_mats[i][j],end=' ')
    print()


