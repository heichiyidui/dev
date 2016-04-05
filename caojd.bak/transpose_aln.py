#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

# mkdir temp
# transpose_aln.py 
# rm -r seq_aln 
# mv temp seq_aln 

import sys

dom_ids=open('index/s.ls').read().split()

for dom_id in dom_ids:
    ifile=open('seq_aln/'+dom_id)
    seqs=[]
    for line in ifile:
        seqs.append(line.strip())
    ifile.close()
    
    ofile=open('temp/'+dom_id,'w')
    ofile.write(str(len(seqs[0]))+' '+str(len(seqs))+'\n')
    for i in range(len(seqs[0])):
        for j in range(len(seqs)):
            ofile.write(seqs[j][i])
        ofile.write('\n')
    ofile.close()
    
