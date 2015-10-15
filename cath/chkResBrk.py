#!/usr/bin/env python3 

import sys
import numpy 

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
    
    C_s  = []
    CA_s = []
    N_s  = []
    
    for line in atom_lines:
        if line[12:16]==' C  ':
            C_s.append( numpy.array([float(line[30:38]),\
                                     float(line[38:46]),\
                                     float(line[46:54])]) )
        if line[12:16]==' CA ':
            CA_s.append(numpy.array([float(line[30:38]),\
                                     float(line[38:46]),\
                                     float(line[46:54])]) )
        if line[12:16]==' N  ':
            N_s.append( numpy.array([float(line[30:38]),\
                                     float(line[38:46]),\
                                     float(line[46:54])]) )
    
    has_residue_backbone_break = False
    backbone_broken_residues=set([])
    for i in range(len(C_s)):
        dis_N_CA = numpy.linalg.norm(CA_s[i]-N_s[i])
        dis_CA_C = numpy.linalg.norm(CA_s[i]-C_s[i])
        
        if dis_N_CA > 2.5 or dis_N_CA < 1.0:
            print(residues[i],"dis_N_CA",dis_N_CA)
            has_residue_backbone_break = True 
            backbone_broken_residues.add(residues[i])
        if dis_CA_C > 2.5 or dis_CA_C < 1.0:
            print(residues[i],"dis_CA_C",dis_CA_C)
            has_residue_backbone_break = True
            backbone_broken_residues.add(residues[i])
    
    if has_residue_backbone_break: 
        print(id,end=' ')
        for res in backbone_broken_residues:
            print(res,end=' ')
        print()

