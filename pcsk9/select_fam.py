#!/usr/bin/env python3

ifile=open('t.dat')
ck_rank = {}
for line in ifile:
    cols = line[:-1].split()
    ck_rank [ cols[1] ] = int(cols[0])
ifile.close()

ifile=open('t_genome.in')
ifile.readline()
to_remove_cks=set([])
for line in ifile:
    cols = line[:-1].split()
    if cols[1] in to_remove_cks:
        continue
    if cols[3] in to_remove_cks:
        continue
    if ck_rank[cols[1]] > ck_rank[cols[3]]:
        to_remove_cks.add(cols[1])
    else:
        to_remove_cks.add(cols[3])
ifile.close()

for ck_id in to_remove_cks:
    print(ck_id)

