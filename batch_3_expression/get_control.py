#!/usr/bin/env python3

s_ids=[]
ifile=open('input_files/batch_info.txt')
ifile.readline()
for line in ifile:
    s_ids.append(line.split()[0])

ifile.close()

ifile=open('input_files/controlprobe.txt.bak')
for i in range(8):
    ifile.readline()
cols=ifile.readline().split('\t')
to_keep_cols=[0,1]
print('controlType'+'\t'+'ProbeID',end='')
for i in range(len(cols)):
    if cols[i].replace('.AVG_Signal','') in s_ids:
        to_keep_cols.append(i)
        print('\t'+cols[i].replace('.AVG_Signal',''),end='')
print()

to_keep_cols=set(to_keep_cols)
for line in ifile:
    cols=line.split('\t')
    
    print(cols[0],end='')
    for i in range(1,len(cols)):
        if i in to_keep_cols:
            print('\t'+cols[i],end='')
    print()
        
