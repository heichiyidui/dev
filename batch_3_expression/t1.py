#!/usr/bin/env python3

ifile=open('t.in')
for line in ifile:
    cols=line.split('\t')
    print(cols[18:-4].count('FALSE'),cols[18:-4].count('TRUE'))
        
    