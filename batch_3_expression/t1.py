#!/usr/bin/env python3 

import numpy

ifile=open('t.in')

predicts=[]
eu_codes=[]
for line in ifile:
    cols=line[:-1].split()
    predicts.append(cols[1])
    eu_codes.append(cols[3])
    
ifile.close()

import random
random.seed(514)
for l in range(100000):
    sum_identical=0
    random.shuffle(predicts)
    
    for i in range(len(eu_codes)):
        for j in range(i+1,len(eu_codes)):
            if eu_codes[i] == eu_codes[j]:
                if predicts[i] == predicts[j]:
                    sum_identical += 1
    print(sum_identical)