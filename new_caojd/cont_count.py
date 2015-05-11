#!/usr/bin/env python3

dom_ids=open('t.ls').read().split()

dom_lens={}
ifile=open('index/cath_s35.dssp')
for line in ifile:
    dom_id=line.strip()[1:]
    line=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    dom_lens[dom_id]=len(line)
ifile.close()

ifile=open('../cath/index/cath_s35.condef')
for line in ifile:
    dom_id=line.strip()[1:]
    line=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    conts=line.split()
    sum_loc_cont=0
    sum_glo_cont=0
    for cont in conts:
        (conti,contj)=cont.split('-')
        if int(contj)-int(conti) < 5:
            sum_loc_cont += 1
        else:
            sum_glo_cont += 1
    print(dom_id,dom_lens[dom_id], sum_loc_cont,sum_glo_cont)
        
