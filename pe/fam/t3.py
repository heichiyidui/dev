#! /usr/bin/env python3
import math

index=open('t.in')
sex={}
age={}
for line in index:
    col=line.split()
    sex[col[0]]=int(col[1])
    age[col[0]]=float(col[2])
index.close()

ifile=open('po.pairs')
for line in ifile:
    col=line.split()
    if age[col[1]] == -9 or age[col[0]] == -9:
        continue
    agediff=abs( age[col[1]] - age[col[0]])
    if agediff < 14:
        print(col[0],sex[col[0]],age[col[0]])
        print(col[1],sex[col[1]],age[col[1]],'\n')
ifile.close()
