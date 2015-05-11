#!/usr/bin/env python3 
##!/home/klinbrc/bin/python3
import sys

package_dir_a='/home/klinbrc/.local/lib/python3.2/site-packages'
sys.path.insert(0, package_dir_a)
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
    
    res_C=[]
    res_N=[]
    for line in atom_lines:
        if line[12:16]==' C  ':
            res_C.append( numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
        if line[12:16]==' N  ':
            res_N.append( numpy.array([float(line[30:38]),\
                                       float(line[38:46]),\
                                       float(line[46:54])]))
    num_brk=0
    for i in range(1,len(residues)):
        if numpy.linalg.norm(res_C[i-1]-res_N[i]) > 2.5:
            print(residues[i],' c-n ',numpy.linalg.norm(res_C[i-1]-res_N[i]))
            num_brk+=1
    if num_brk>0:
        print('file:',id,len(residues),num_brk)

