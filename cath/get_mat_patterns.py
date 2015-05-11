#!/usr/bin/env python3 

PATTERN_WID = 5

ofile_name='t00.dat'
ofile=open(ofile_name,'w')
for i in range(PATTERN_WID):
    for j in range(PATTERN_WID):
        ofile.write('0 ')
    ofile.write('\n')
ofile.close()

ifile=open('index/pat_67.txt')
ifile.readline()
clu_num=0
for line in ifile:
    clu_num +=1
    ofile_name = 't%02d.dat' % clu_num
    cols=line.split()
    
    ofile=open(ofile_name,'w')
    col_num = 0
    for i in range(PATTERN_WID):
        for j in range(PATTERN_WID):
            ofile.write(cols[col_num]+' ')
            col_num += 1
        ofile.write('\n')
    ofile.close()
ifile.close()
