#!/usr/bin/env python3

ifile=open('t.in')
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for id in domain_ids:
    ifile=open('dompdb/'+id)
    atom_ids=set()       # resName, chainId, resNumber, insertCode and atomName
    atom_min_altCode={}  # smallest altCode for this atom 
    for line in ifile:
        atom_id=line[17:20]+line[21]+line[22:26]+line[26]+line[12:16]
        if atom_id not in atom_ids:
            atom_ids.add(atom_id)
            atom_min_altCode[atom_id]=line[16]
        else:
            if line[16] < atom_min_altCode[atom_id]:
                atom_min_altCode[atom_id]=line[16]
    ifile.seek(0)
    
    ofile=open('pdb2/'+id,'w')
    for line in ifile:
        atom_id=line[17:20]+line[21]+line[22:26]+line[26]+line[12:16]
        if line[16] == atom_min_altCode[atom_id]:
            ofile.write(line)
    ofile.close()
    ifile.close()


#0         1         2         3         4         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#ATOM     86  CG  ARG    11      -2.455   1.706  24.211  1.00 17.72      1AAK 146
#ATOM      1  N   GLY A   0       3.756  -1.147  19.133  1.00 45.33
