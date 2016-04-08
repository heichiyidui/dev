#!/usr/bin/env python3
import sys
batch_id = sys.argv[1]
# batch_id = 'b01'

ifile=open(batch_id+'/posterior.txt')
ifile.readline()
print('id x1 vx1 y1 vy1 cov1 x2 vx2 y2 vy2 cov2 x3 vx3 y3 vy3 cov3')

for line in ifile:
    cols = line[:-1].split('\t')
    snp_id = cols[0]
    bb_cols = cols[1].split(',')
    ab_cols = cols[2].split(',')
    aa_cols = cols[3].split(',')

    print(snp_id,end=' ')
    print(bb_cols[0],bb_cols[1],bb_cols[4],bb_cols[5],bb_cols[6],end=' ')
    print(ab_cols[0],ab_cols[1],ab_cols[4],ab_cols[5],ab_cols[6],end=' ')
    print(aa_cols[0],aa_cols[1],aa_cols[4],aa_cols[5],aa_cols[6],end='\n')
ifile.close()