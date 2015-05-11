#!/usr/bin/env python3 
##!/home/klinbrc/bin/python3
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
    
    hasC= [False]*len(residues)
    hasCA=[False]*len(residues)
    hasN= [False]*len(residues)
    
    for line in atom_lines:
        if line[12:16]==' C  ':
            hasC[residues.index(line[21:27])]=True
        if line[12:16]==' CA ':
            hasCA[residues.index(line[21:27])]=True
        if line[12:16]==' N  ':
            hasN[residues.index(line[21:27])]=True
    
    brk_res_set=set([])
    for i in range(len(residues)):
        if hasC[i] == False or hasCA[i] == False or hasN[i] == False:
            brk_res_set.add(residues[i])
    
    print(id,len(residues),\
        hasC.count(False),hasCA.count(False),hasN.count(False),len(brk_res_set))
    
    if len(brk_res_set)>0:
        ofile=open('pdb2/'+id,'w')
        for line in atom_lines:
            if line[21:27] not in brk_res_set:
                ofile.write(line+'\n')
        ofile.close()


