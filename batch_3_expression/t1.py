#!/usr/bin/env python3 
probe_set=[]

ifile=open('50_probes.ls')
for line in ifile :
    probe_set.append(line[:-1])
ifile.close()

#[11]
ifile=open('/home/kuang/dev/batch_3_expression/ba3/BA3_lumi_processing_t_out/'+\
    'BA3.eset_bg_log2_rsn_adj.exprs_matrix.txt')
line1=ifile.readline()[:-1]
cols=line1.split('\t')
print('PROBE_ID ' + ' '.join(cols[30:]))

for line in ifile:
    cols=line[:-1].split('\t')
    if cols[11] not in probe_set:
        continue
    print(cols[11],end=' ')
    print(' '.join(cols[30:]))
ifile.close()