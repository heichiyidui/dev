#!/usr/bin/env python3 
import sys 
import numpy 

##!/home/klinbrc/bin/python3
#import sys
#package_dir_a='/home/klinbrc/.local/lib/python3.2/site-packages'
#sys.path.insert(0, package_dir_a)
#import numpy

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for id in domain_ids:
    ifile=open('dompdb/'+id)
    atom_lines=[]
    atom_xyzs=[]
    for line in ifile:
        xyz=numpy.array([float(line[30:38]),\
                         float(line[38:46]),\
                         float(line[46:54])])
        atom_xyzs.append(xyz)
        atom_lines.append(line)
    ifile.close()
    
    has_close_atoms=False;
    domain_size=len(atom_xyzs)
    for i in range(domain_size):
        for j in range(i+1,domain_size):
            dis_ij = numpy.linalg.norm(atom_xyzs[i]-atom_xyzs[j])
            if dis_ij < 0.5 :
                print(atom_lines[i]+atom_lines[j],end='')
                print('distance: ','%6.3f' % dis_ij)
                has_close_atoms=True
    
    if has_close_atoms:
        print('file:',id,'\n')
