#!/usr/bin/env python3

AA_NAME_TO_CODE={'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',\
                 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',\
                 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',\
                 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

dom_ids=open('t.ls').read().split()

seqs={}
ifile=open('index/cath_s35.seq')
for line in ifile:
    dom_id=line.strip()[1:]
    seqs[dom_id]=ifile.readline().strip()
ifile.close()

for dom_id in dom_ids:
    seq=seqs[dom_id]
    
    acc=[]
    ss=[]
    ifile=open('stride_ss/'+dom_id)
    seq_index=0
    for line in ifile:
        if not line.startswith('ASG') :
            continue
        if AA_NAME_TO_CODE[line[5:8]] != seq[seq_index]:
            print(dom_id,seq_index,line.strip())
        seq_index += 1
        ss.append(line[24])
        acc.append(line[61:69].strip())
    ifile.close()
    
    print('>'+dom_id)
    # print(''.join(ss))
    # print(' '.join(acc))
    
    
    

