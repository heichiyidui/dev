#!/usr/bin/env python3

cath_class = {}
ifile=open('../cath/index/CathDomainList.S35')
for line in ifile:
    cols = line[:-1].split()
    dom_class = '_'.join(cols[1:5])
    dom_id    = cols[0]
    cath_class[dom_id] = dom_class
ifile.close()

collected_classes=set([])
dom_reprs=[]
ifile=open('t.in')
for line in ifile:
    cols = line[:-1].split()
    aln_size,dom_id = int(cols[0]) , cols[1]
    if cath_class[dom_id] in collected_classes:
        continue
    collected_classes.add(cath_class[dom_id])
    dom_reprs.append(dom_id)
ifile.close()

for dom_id in dom_reprs:
    print(dom_id)

