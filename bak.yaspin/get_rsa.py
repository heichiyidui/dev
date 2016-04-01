#!/usr/bin/env python3

dom_ids=open('index/s.ls').read().split()

cseqs={}
ifile=open('../cath/index/cath_s35.cseq')
for line in ifile:
    dom_id=line.strip()[1:]
    cseq=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    cseqs[dom_id]=cseq 

ifile.close()

# max residue acc Values from Sander & Rost, (1994), Proteins, 20:216-226
MAX_ACC={}
MAX_ACC['A']=106.0      
MAX_ACC['R']=248.0      
MAX_ACC['N']=157.0      
MAX_ACC['D']=163.0      
MAX_ACC['C']=135.0      
MAX_ACC['Q']=198.0      
MAX_ACC['E']=194.0      
MAX_ACC['G']= 84.0      
MAX_ACC['H']=184.0      
MAX_ACC['I']=169.0      
MAX_ACC['L']=164.0      
MAX_ACC['K']=205.0      
MAX_ACC['M']=188.0      
MAX_ACC['F']=197.0      
MAX_ACC['P']=136.0      
MAX_ACC['S']=130.0      
MAX_ACC['T']=142.0      
MAX_ACC['W']=227.0      
MAX_ACC['Y']=222.0      
MAX_ACC['V']=142.0    

ifile=open('../cath/index/cath_s35_stride.acc')
for line in ifile:
    dom_id = line.strip()[1:]
    accs=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    accs=accs.split()

    for i in range(len(accs)):
        if cseqs[dom_id][i] not in MAX_ACC.keys():
            continue
        print(float(accs[i])/MAX_ACC[cseqs[dom_id][i]])

ifile.close()
