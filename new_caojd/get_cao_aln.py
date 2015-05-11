#!/home/klinbrc/bin/python3 

import sys 

MAX_ALN_DEPTH  = 2000

AA_CODE_TO_INT = {'A': 0 ,'R': 1 ,'N': 2 ,'D': 3 ,'C': 4 ,\
                  'Q': 5 ,'E': 6 ,'G': 7 ,'H': 8 ,'I': 9 ,\
                  'L': 10,'K': 11,'M': 12,'F': 13,'P': 14,\
                  'S': 15,'T': 16,'W': 17,'Y': 18,'V': 19,};

dom_ids=open(sys.argv[1]).read().split()

dom_depths={}
ifile=open('index/cath_s35.aln_depth')
for line in ifile:
    cols=line.split()
    if cols[0] not in dom_ids:
        continue
    dom_depths[cols[0]]=int(cols[1])
ifile.close()

# read domain contact definitions
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

import random
random.seed(514)

# read and write alignments    
for dom_id in dom_ids:
    ifile=open('compact_paln/'+dom_id)
    ifile.readline()
    dom_line=ifile.readline().strip()
    dom_seq=[]
    for c in dom_line:
        dom_seq.append(AA_CODE_TO_INT[c]) 
    dom_len=len(dom_seq)
    
    ofile=open('cao_aln/'+dom_id,'w')
    for line in ifile:
        subjc_seq=ifile.readline().strip()
        
        if dom_depths[dom_id] > MAX_ALN_DEPTH:
            if random.random() > MAX_ALN_DEPTH/dom_depths[dom_id]:
                continue
        
        aln=[0]*160000
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
            cont_dom   = dom_seq[c_i]*20 + dom_seq[c_j]
            cont_subjc = AA_CODE_TO_INT[subjc_seq[c_i]]*20 \
                        +AA_CODE_TO_INT[subjc_seq[c_j]]
            
            aln[cont_dom*400 + cont_subjc] += 1
            aln[cont_subjc*400 + cont_dom] += 1
            
        ofile.write(' '.join(str(x) for x in aln)+'\n')
    ofile.close()
    ifile.close()
         

