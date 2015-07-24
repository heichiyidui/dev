#!/usr/bin/env python3 

import numpy

ifile=open('t.in')

predicts=[]
eu_codes=[]
for line in ifile:
    cols=line[:-1].split()
    if cols[3] != "MCI":
        continue
    predicts.append(cols[1])
    eu_codes.append(cols[2])
ifile.close()

to_be_compaired_ij = set()
for i in range(len(eu_codes)):
    for j in range(i+1,len(eu_codes)):
        if eu_codes[i] == eu_codes[j]:
            to_be_compaired_ij.add((i,j))

# print(len(to_be_compaired_ij))
# sum_identical=0
# for ij in to_be_compaired_ij:
#     if predicts[ij[0]] == predicts[ij[1]]:
#         sum_identical += 1
# print(sum_identical)
# #############213 identical out of 348 comparisons

import random
random.seed(514)

for l in range(100000):
    sum_identical = 0 
    random.shuffle(predicts)
    for ij in to_be_compaired_ij:
        if predicts[ij[0]] == predicts[ij[1]]:
            sum_identical += 1
    print(sum_identical)


