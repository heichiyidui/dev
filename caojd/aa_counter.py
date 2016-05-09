#!/usr/bin/env python3
import collections

aa_counter = collections.Counter()

dom_seq={}
ifile=open('../cath/index/cath_s35.seq')
for line in ifile:
    dom_id = line[1:-1]
    seq = ifile.readline()[:-1]
    dom_seq[dom_id] = seq
ifile.close()

dom_ls = open('index/dom.ls').read().split()
for dom_id in dom_ls:
    aa_counter.update(dom_seq[dom_id])

aa_sum = sum(aa_counter.values())
for c in aa_counter.keys():
    print(c,aa_counter[c],'{:8.5f}'.format(aa_counter[c]/aa_sum) )
