#! /usr/bin/env python3

ifile=open('t.ls')
ind=[];
for line in ifile:
    ind.append(line.strip())
ifile.close()

ifile=open('pe8.fam')
for line in ifile:
    col=line.split()
    if col[1] not in ind:
        continue;
    print(line[:-1])
ifile.close()
