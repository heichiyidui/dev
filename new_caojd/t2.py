#!/usr/bin/env python3 

dom_ids=open('t.ls').read().split()

ss7_lines={}
ifile=open('index/cath_s35.stride.ss7')
for line in ifile:
    dom_id=line[1:-1]
    ss_line=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    ss7_lines[dom_id]=ss_line 
ifile.close()

res_count=0
ifile=open('index/cath_s35.dssp')
for line in ifile:
    dom_id=line[1:-1]
    ss_line=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    if len(ss_line) != len(ss7_lines[dom_id]):
        print(dom_id)
    res_count+=len(ss_line)
ifile.close()

print(res_count)