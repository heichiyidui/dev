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
    
    residues=[]
    residue_set=set([])
    for line in atom_lines:
        if line[21:27] not in residue_set:
            residues.append(line[21:27])
            residue_set.add(line[21:27])
    
    res_cnca_alt={} # altLoc codes of C,CA,N atoms
    for res in residues:
        res_cnca_alt[res]=[]
    
    cnca_set=set([' C  ',' CA ',' N  '])
    for line in atom_lines:
        if line[12:16] not in cnca_set:
            continue
        if line[16] == ' ':
            continue 
        res_cnca_alt[line[21:27]].append(line[16])
    
    has_multi_alt_cnca=False
    for res in residues:
        if len(set(res_cnca_alt[res])) > 1:
            print(res,res_cnca_alt[res])
            has_multi_alt_cnca=True
    
    if has_multi_alt_cnca:
        print('file:',id)
        
            
            
