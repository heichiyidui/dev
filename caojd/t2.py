#!/usr/bin/env python3

dom_ids=open('index/s.ls').read().split()

for dom_id in dom_ids:
    ifile=open('seq_aln/'+dom_id)
    cols=ifile.readline().split()
    dom_len=int(cols[0])
    aln_depth=int(cols[1])
    
    for i in range(dom_len):
        line=ifile.readline().strip()
        C0=line[0]
        print(line.count(C0)/aln_depth)
    ifile.close()
    