#!/usr/bin/env python3

ifile=open('index/CathDomainList.S35')
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for id in domain_ids:
    ifile=open('dompdb/'+id)
    atom_ids=set()    # resName, chainId, resNumber, insertCode and atomName
    multiloc_res=set()# resName, chainId, resNumber, insertCode
    for line in ifile:
        atom_id=line[17:20]+line[21]+line[22:26]+line[26]+line[12:16]
        if atom_id in atom_ids:
            multiloc_res.add(line[17:20]+line[21]+line[22:26]+line[26])
        atom_ids.add(atom_id)
    ifile.close()
    if len(multiloc_res) > 0:
        print(id,len(multiloc_res))


#0         1         2         3         4         5         6         7
#012345678901234567890123456789012345678901234567890123456789012345678901
#ATOM     86  CG  ARG    11      -2.455   1.706  24.211  1.00 17.72      
#ATOM      1  N   GLY A   0       3.756  -1.147  19.133  1.00 45.33
