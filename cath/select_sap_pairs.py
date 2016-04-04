#!/usr/bin/env python3

ifile=open('index/CathDomainList.S35')
cath_class={}
for line in ifile:
    cols = line[:-1].split()
    dom_id = cols[0]
    cath   = '_'.join(cols[1:5])
    if cath not in cath_class.keys():
        cath_class[cath]=[dom_id]
    else:
        cath_class[cath].append(dom_id)
ifile.close()
print(len(cath_class))

keyset=set(cath_class.keys())
for cath in keyset:
    if len(cath_class[cath]) < 2:
        del cath_class[cath]

print(len(cath_class))

import random
random.seed(514)

pair_set=set()
for i in range(40):
    for cath in cath_class.keys():
        random.shuffle(cath_class[cath])
        pair=sorted(cath_class[cath][0:2])
        pair_set.add(' '.join(pair))

# for p in pair_set:
#     id1,id2=p.split()
#     print('~/bin/sap','dompdb/'+id1,'dompdb/'+id2,'> sap_aln/'+id1+'_'+id2)

# 1427 pairs in the first iteration
#  836 in the second
#  762 in the third ...

