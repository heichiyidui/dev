#!/usr/bin/env python3 

ifile=open('t.in')
for line in ifile:
    cols=line[:-1].split()
    if len(cols) == 2:
        print(cols[0],cols[1])
    else: 
        print(cols[0],cols[2])
        print(cols[1],cols[2])

ifile.close()