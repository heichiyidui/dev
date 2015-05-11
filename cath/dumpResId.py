#!/usr/bin/env python3 
import sys
import pickle

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
        res_id_map[new_res]=residues[i]
    
    ofile=open('res_num_map/'+id,'wb')
    pickle.dump(res_id_map,ofile)
    ofile.close()

