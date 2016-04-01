#!/usr/bin/env python3

written_cathh_ids=[]
ifile=open('index/domain.ls')
for line in ifile:
    cols=line.split()
    id=cols[0]
    cathh_id=cols[1]+'.'+cols[2]+'.'+cols[3]+'.'+cols[4]
    if cathh_id not in written_cathh_ids:
        print(line[:-1])
        written_cathh_ids.append(cathh_id)
        