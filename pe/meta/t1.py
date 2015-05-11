#! /usr/bin/env python3

ifile=open('t.in')
ids=[]
for line in ifile:
    if line.split()[1] == 'WAFSS':
#    if line.split()[1] == 'PE_PAMP' or line.split()[1] == 'Valdecilla' :
        ids.append(line.split()[0])
ifile.close()

ifile=open('../assoc/pe17.fam')
for line in ifile:
    if line.split()[1] in ids:
        print(line[:-1])
ifile.close()
