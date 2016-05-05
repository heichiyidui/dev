#!/usr/bin/env python3

import sys

families = []

ifile=open(sys.argv[1])
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    if float(cols[9]) < 0.125:
        continue
    id1=cols[1]
    id2=cols[3]

    is_pair_found = []
    for i in range(len(families)):
        if id1 in families[i] or id2 in families[i]:
            is_pair_found.append(i)

    if is_pair_found == []:
        families.append(set([id1,id2]))
        continue
    if len(is_pair_found) == 1:
        families[is_pair_found[0]].add(id1)
        families[is_pair_found[0]].add(id2)
        continue

    if len(is_pair_found) > 1: # merge sets
        new_id_set = set([])
        for i in is_pair_found:
            for id1 in families[i]:
                new_id_set.add(id1)
        new_families = []
        for i in range(len(families)):
            if i not in is_pair_found:
                new_families.append(families[i])
        new_families.append(new_id_set)
        families = new_families

ifile.close()

for i in range(len(families)):
    print(len(families[i]))

