#!/usr/bin/env python3 

import numpy 

ifile=open('input_files/ba3_50_probes.txt')
print(ifile.readline()[:-1])
for line in ifile:
    cols=line[:-1].split()
    print(cols[0],end=' ')
    
    expres_prof=[]
    for i in range(1,len(cols)):
        expres_prof.append(float(cols[i]))
    
    mean = numpy.average(expres_prof)
    std  = numpy.std(expres_prof)
    
    out_cols=[]
    for i in range(len(expres_prof)):
        expres_prof[i] -= mean
        expres_prof[i] /= std 
        out_cols.append('{:6.3f}'.format(expres_prof[i]))
    print(' '.join(out_cols))
    
ifile.close()