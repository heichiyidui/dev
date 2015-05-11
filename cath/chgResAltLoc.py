#!/usr/bin/env python3 
import sys

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for id in domain_ids:
    # reading the atom lines
    ifile=open('pdb/'+id)
    atom_lines=ifile.read().split('\n')[:-1]
    ifile.close()
    
    # reading residue set
    residues=[] 
    residue_set=set([])
    for line in atom_lines:
        if line[21:27] not in residue_set:
            residues.append(line[21:27])
            residue_set.add(line[21:27])
    
    res_alt_codes={} 
    for res in residues:
        res_alt_codes[res]=[]
    
    for line in atom_lines:
        if line[16] == ' ':
            continue 
        res_alt_codes[line[21:27]].append(line[16])
    
    # finding the residues with multiple altloc codes
    multi_alt_loc_res=set([])
    max_alt_loc_code={}
    for res in residues:
        code_list = res_alt_codes[res]
        if len(set(code_list)) < 2:
            continue
        multi_alt_loc_res.add(res)
        max_code = max(set(code_list),key=code_list.count)
        max_alt_loc_code[res] = max_code
    # now we know which residues to change and which codes to use
    
    ifile=open('rscb/'+id[:4]+'.pdb')
    pdb_lines=ifile.read().split('\n')[:-1]
    ifile.close()
    
    ofile=open('pdb2/'+id,'w')
    for line in atom_lines:
        res=line[21:27]
        if res not in multi_alt_loc_res:
            ofile.write(line+'\n')
            continue
        if line[16] ==' ' or line[16] == max_alt_loc_code[res]:
            ofile.write(line+'\n')
            continue
        
        is_line_replaced=False
        for p_line in pdb_lines:
            if p_line[21:27] != res:
                continue
            if p_line[12:16] != line[12:16]: # atom name 
                continue
            if p_line[16] == max_alt_loc_code[res]:
                ofile.write(p_line+'\n')
                is_line_replaced=True
                break
        
        if not is_line_replaced:
            if line[12:16].startswith(' H'):
                continue
            if line[12:16].startswith('H'):
                continue
            
            ofile.write(line+'\n')
            print(id,res,'\'',max_alt_loc_code[res],'\'',line)
    
    ofile.close()
    print('')    

