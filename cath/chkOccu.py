#!/usr/bin/env python3

ifile=open('index/CathDomainList.S35')
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for id in domain_ids:
    ifile=open('dompdb/'+id)
    for line in ifile:
        if float(line[54:60]) < 0:
            print(id,line[:-1])
    ifile.close()

