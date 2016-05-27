#!/usr/bin/env python3

ifile=open('t.dat')
ck_links = {}
for line in ifile:
    cols = line[:-1].split()
    ck_links [ cols[1] ] = int(cols[0])
ifile.close()

ifile=open('t.in')
ifile.readline() # the header line
to_remove_cks=set([])
for line in ifile:
    cols = line[:-1].split()
    if cols[0] in to_remove_cks:
        continue
    if cols[1] in to_remove_cks:
        continue
    if ck_links[cols[0]] > ck_links[cols[1]]:
        to_remove_cks.add(cols[0])
    else:
        to_remove_cks.add(cols[1])
ifile.close()

for ck_id in to_remove_cks:
    print(ck_id)

