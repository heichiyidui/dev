#!/usr/bin/env python3 

##!/home/klinbrc/bin/python3
# package_dir_a='/home/klinbrc/.local/lib/python3.2/site-packages'
# sys.path.insert(0, package_dir_a)

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
    
    atom_xyzs=[]
    for line in atom_lines:
        xyz=numpy.array([float(line[30:38]),\
                         float(line[38:46]),\
                         float(line[46:54])])
        atom_xyzs.append(xyz)
    
    has_brk_res=False;
    domain_size=len(atom_xyzs)
    for i in range(domain_size):
        for j in range(i+1,domain_size):
            if atom_lines[i][21:27] != atom_lines[j][21:27]:
                continue
            dis_ij = numpy.linalg.norm(atom_xyzs[i]-atom_xyzs[j])
            if dis_ij > 11 :
                print(atom_lines[i]+'\n'+atom_lines[j])
                print('distant:','%6.3f' % dis_ij)
                has_brk_res=True
    
    if has_brk_res:
        print('file:',id,'\n')
