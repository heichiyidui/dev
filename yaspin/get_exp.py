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

import math 
pi1 = 0.3378852 ; k1 = 0.716621 ; l1 = 0.1162887 
pi2 = 0.1152352 ; k2 = 1.904627 ; l2 = 0.2615063 
pi3 = 0.5468796 ; k3 = 2.664756 ; l3 = 0.6059084 

ifile=open('../cath/index/cath_s35_stride.acc')
for line in ifile:
    dom_id = line.strip()[1:]
    accs=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    accs=accs.split()
    
    rsas=[] # relative solvent access
    for i in range(len(accs)):
        if accs[i]=='NA':
            rsas.append('NA')
            continue
        c = cseqs[dom_id][i]
        rsas.append(float(accs[i])/MAX_ACC[c])

    exps=[] # residue exposure
    for i in range(len(rsas)):
        if rsas[i] == 'NA':
            exps.append('NA')
            continue
        x=rsas[i]
        if x==0:
            exps.append(0.0)
            continue
        if x>0.75:
            exps.append(1.0)
            continue

        p1 = pi1 * (k1/l1) * ((x/l1) ** (k1-1)) *math.exp(- (x/l1)**k1 )
        p2 = pi2 * (k2/l2) * ((x/l2) ** (k2-1)) *math.exp(- (x/l2)**k2 )
        p3 = pi3 * (k3/l3) * ((x/l3) ** (k3-1)) *math.exp(- (x/l3)**k3 )
        exps.append(p3/(p1+p2+p3))

    print('>'+dom_id)
    for exp in exps:
        print(exp,end=' ')
    print()

ifile.close()
