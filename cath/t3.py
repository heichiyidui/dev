#!/usr/bin/env python3 
##!/home/klinbrc/bin/python3

PATTERN_WID = 25

# read clustering centers from the last clustering
clus=[int(clu) for clu in open('t.clu').read().split()]

CLU_NUM = max(clus)+1

# sum sub-matrices to patterns

patterns = []
for i in range(CLU_NUM): patterns.append([0]*PATTERN_WID)
pattern_sums=[0]*CLU_NUM

lines=open('t.in').read().split('\n')
while '' in lines: lines.remove('')
for i in range(len(clus)):
    line=lines[i]
    sub_mat = [float(col) for col in line.split()]

    for j in range(PATTERN_WID):
        patterns[clus[i]][j] += sub_mat[j]
    pattern_sums[clus[i]] += 1

for i in range(CLU_NUM):
    for j in range(PATTERN_WID):
        print(patterns[i][j] / pattern_sums[i],end=' ')
    print()

