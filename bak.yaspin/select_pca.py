#!/usr/bin/env python3
import math
import random

WING_SIZE=1

# Values from Sander & Rost, (1994), Proteins, 20:216-226
# Used for relative accessibility
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

# we wants to ignore the few 'I' and select all B and keep the other states
# about the same number with B
select_ss={}
select_ss['B']=1.0                
select_ss['C']=0.04876314620056288
select_ss['E']=0.04944576286460993
select_ss['G']=0.280289484887186  
select_ss['H']=0.02978646398841838
select_ss['I']=1.0 # should be 39.7
select_ss['S']=0.12350403301444382
select_ss['T']=0.09450265537534089


random.seed(514)

ids=[]
ifile=open('index/h.ls')
for line in ifile:
    ids.append(line.split()[0])
ifile.close()

for id in ids:
    ifile=open('in/'+id)
    lines = ifile.read().split('\n')[:-1]
    ifile.close()
    
    for i in range(WING_SIZE,len(lines)-WING_SIZE):
        
        ss=lines[i][4]
        if ss=='-':
            continue
        
        t_f=random.random()
        if t_f > select_ss[ss]:
            continue
        
        print(ss+' ',end='')
        
        aa=lines[i][2]
        rac=float(lines[i][6:10])/MAX_ACC[aa]
        print(rac,' ',end='')
#        if random.random() < 0.1 :
#            print(rac)
        
        
#        for j in range(-WING_SIZE,WING_SIZE+1):
#            line=lines[i+j]
#            for k in range(20):
#                print(line[11+k*3:14+k*3]+' ',end='')
#        print()
#            
#    break;
    
