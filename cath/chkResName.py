#!/usr/bin/env python3 
import sys

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for id in domain_ids:
    ifile=open('dompdb/'+id)
    atom_lines=ifile.read().split('\n')[:-1]
    ifile.close()
    
    res_atom_list=[]
    for line in atom_lines:
        res_atom_list.append((line[21:27],line))
    
    res_atom_list.sort()
    
    hasMultiResName=False;
    for i in range(1,len(atom_lines)):
        if res_atom_list[i-1][0] != res_atom_list[i][0]:
            continue
        atom1=res_atom_list[i-1][1]
        atom2=res_atom_list[i]  [1]
        if atom1[17:20] != atom2[17:20]:
            print(atom1+'\n'+atom2)
            hasMultiResName=True
    
    if hasMultiResName:
        print('file:',id)
    
    

