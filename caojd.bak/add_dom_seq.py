#!/usr/bin/env python3
import sys
dom_ids=open(sys.argv[1]).read().split()

ifile=open('index/cath_s35.seq')
dom_seq={}
for line in ifile:
    s_line=ifile.readline().strip()
    dom_id=line.strip()[1:]
    if dom_id not in dom_ids:
        continue
    dom_seq[dom_id]=s_line
ifile.close()

for dom_id in dom_ids:
    ifile=open('bseqs/'+dom_id)
    ofile=open('temp/'+dom_id,'w')
    ofile.write('>'+dom_id+'\n')
    ofile.write(dom_seq[dom_id]+'\n')
    for line in ifile:
        if line.startswith('>'):
            cols=line.split()
            gi_id=cols[0].split('|')[1]
            gi_range=cols[0].split(':')[-1]
            ofile.write('>gi|'+gi_id+'|:'+gi_range+'\n')
        else:
            ofile.write(line)
    
    ifile.close()            
    ofile.close()
            