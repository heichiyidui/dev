#!/usr/bin/env python3

AA_NAME_TO_CODE={'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',\
                 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',\
                 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',\
                 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

dom_ids = open('t.ls').read().split()

for dom_id in dom_ids:
    ifile=open('dompdb/'+dom_id)
    atom_lines=ifile.read().split('\n')[:-1]
    ifile.close()

    seq=''
    for line in atom_lines:
        if line[12:16]==' CA ':
            seq+=AA_NAME_TO_CODE[line[17:20]]
    print(">"+dom_id+'\n'+seq)

