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
    
    # read the residue set
    residues=[] 
    for line in atom_lines:
        if line[21:27] not in residues:
            residues.append(line[21:27])
    
    res_id_map={}
    for i in range(len(residues)):
        new_res=' %4d '%(i+1)
        res_id_map[residues[i]]=new_res
    
    ofile=open('pdb2/'+id,'w')
    for line in atom_lines:
        ofile.write(line[:16]+' '+line[17:21]\
                   +res_id_map[line[21:27]]+line[27:]+'\n')
    ofile.close()
    
