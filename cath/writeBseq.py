#!/usr/bin/env python3 
import sys
import numpy

AA_NAME_TO_CODE={'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',\
                 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',\
                 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',\
                 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}


ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for id in domain_ids:
    ifile=open('dompdb/'+id)
    atom_lines=ifile.read().split('\n')[:-1]
    ifile.close()
    
    seq=''
    for line in atom_lines:
        if line[12:16]==' CA ':
            seq+=AA_NAME_TO_CODE[line[17:20]]
    
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
    out_seq =  ''
    out_seq += seq[0]
    for i in range(1,len(seq)):
        if numpy.linalg.norm(res_C[i-1]-res_N[i]) > 2.5:
            out_seq += 'x'
        out_seq+=seq[i]
    
    ofile=open('bseq/'+id,'w')
    ofile.write('>'+id+'\n')
    ofile.write(out_seq)
    ofile.close()
    
