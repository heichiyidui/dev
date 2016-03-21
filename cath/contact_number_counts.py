#!/usr/bin/env python3

ifile = open('index/cath_s35.seq')
seqs={}
for line in ifile:
    dom_id = line[1:-1]
    seq = ifile.readline()[:-1]
    seqs[dom_id] = seq
ifile.close()

ifile=open('index/cath_s35.condef')
for line in ifile:
    dom_id = line[1:-1]
    line = ifile.readline()
    contacts = line[:-1].split()
    res_conts=[0]*len(seqs[dom_id])

    for contact in contacts:
        (a,b) = contact.split('-')
        (cont_i,cont_j) = (int(a),int(b))
        res_conts[cont_i] += 1
        res_conts[cont_j] += 1

    for i in range(len(seqs[dom_id])):
        print(seqs[dom_id][i],res_conts[i])


ifile.close()