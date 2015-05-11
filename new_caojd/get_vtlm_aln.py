#!/home/klinbrc/bin/python3 

MAX_ALN_DEPTH  = 2000

AA_CODE_TO_INT = {'A': 0 ,'R': 1 ,'N': 2 ,'D': 3 ,'C': 4 ,\
                  'Q': 5 ,'E': 6 ,'G': 7 ,'H': 8 ,'I': 9 ,\
                  'L': 10,'K': 11,'M': 12,'F': 13,'P': 14,\
                  'S': 15,'T': 16,'W': 17,'Y': 18,'V': 19,};

import sys                   
dom_ids=open(sys.argv[1]).read().split()
dom_depths={}
ifile=open('index/cath_s35.aln_depth')
for line in ifile:
    cols=line.split()
    if cols[0] not in dom_ids:
        continue
    dom_depths[cols[0]]=int(cols[1])
ifile.close()

import random
random.seed(514)

for dom_id in dom_ids:
    ifile=open('compact_paln/'+dom_id)
    ofile=open('vtml_aln/'+dom_id,'w')
    ifile.readline()
    dom_seq=ifile.readline().strip()
    dom_len=len(dom_seq)
    
    for line in ifile:
        subjc_seq=ifile.readline().strip()
        
        if dom_depths[dom_id] > MAX_ALN_DEPTH:
            if random.random() > MAX_ALN_DEPTH/dom_depths[dom_id]:
                continue
        
        aln=[0]*400
        for i in range(dom_len):
            if subjc_seq[i] == '-':
                continue
            dom_i = AA_CODE_TO_INT[dom_seq[i]]
            sub_i = AA_CODE_TO_INT[subjc_seq[i]]
            
            aln[dom_i*20+sub_i] += 1
            aln[sub_i*20+dom_i] += 1
        ofile.write(' '.join(str(x) for x in aln)+'\n')
    ifile.close()
    ofile.close()     

        