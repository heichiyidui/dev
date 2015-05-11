#!/usr/bin/env python3 

ifile=open('t.in')
for line in ifile:
    cols=line.split()
    x=float(cols[0])
    y=float(cols[1])
    z=float(cols[2])
    group=cols[3]
    
    x=x*10
    y=y*10
    z=z*10
    
    label='  1'
    if group=='1':
        label='100'
    if group=='2':
        label='20O'
    
    print('ATOM   4760  CA  ILE   '+label+'    ',end='')
    print('{:8.3f}{:8.3f}{:8.3f}'.format(x,y,z))
ifile.close()