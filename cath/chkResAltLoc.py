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
    
    residue_set=set([])
    for line in atom_lines:
        residue_set.add(line[21:27])
    
    res_alt_codes={} 
    for res in residue_set:
        res_alt_codes[res]=[]
    
    for line in atom_lines:
        if line[16] == ' ':
            continue 
        res_alt_codes[line[21:27]].append(line[16])
    
    alt_loc_res=set()
    for res in residue_set:
        if len(set(res_alt_codes[res])) > 1:
            alt_loc_res.add(res)
            print(res,res_alt_codes[res])
            print('grep "ATOM............[B]" rcsb/'+id[0:4]+'.pdb',\
                ' | grep "'+res+'"' )
            print('nedit -do \'find(\"'+res+'\")\' pdb2/'+id,"&\n")
    
    if len(alt_loc_res) > 0:
        print('nedit rcsb/'+id[0:4]+'.pdb',"& \n")

            
            
