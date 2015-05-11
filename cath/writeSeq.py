#!/usr/bin/env python3 
import sys

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

AA_NAME_TO_CODE={'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',\
                 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',\
                 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',\
                 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

SEQ_LINE_WIDTH=80

for id in domain_ids:
    ifile=open('dompdb/'+id)
    atom_lines=ifile.read().split('\n')[:-1]
    ifile.close()
    
    seq=''
    for line in atom_lines:
        if line[12:16]==' CA ':
            seq+=AA_NAME_TO_CODE[line[17:20]]
    
    ofile=open('seq/'+id,'w')
    ofile.write('>'+id+'\n')
    for i in range(0,len(seq),SEQ_LINE_WIDTH):
        ofile.write(seq[i:i+SEQ_LINE_WIDTH]+'\n')
    ofile.close()
    
