#!/usr/bin/env python3 
##!/home/klinbrc/bin/python3

PATTEN_WID = 25

# read cluster centers from the last clustering
clus=[]
ifile=open('t.pat')
for line in ifile:
    clu=[float(col) for col in line.split()]
    clus.append(clu)

ifile.close()

CLU_NUM = len(clus)

new_clus=[]
for i in range(CLU_NUM): new_clus.append([0]*PATTEN_WID)
new_clu_sums=[0]*CLU_NUM

# get new centers
import math
sum_dis = 0

lines=open('fxf.in').read().split('\n')
while '' in lines: lines.remove('')
for line in lines:
    sub_mat=[float(col) for col in line.split()]

    dis=[1000]*CLU_NUM
    for i in range(CLU_NUM):
        t_dis=0
        for j in range(PATTEN_WID):
            t_dis += math.fabs(sub_mat[j]-clus[i][j]) # manhattan distance
        dis[i]=t_dis

    clu_index=dis.index(min(dis))
    sum_dis += min(dis)

    new_clu_sums[clu_index] += 1
    for i in range(PATTEN_WID):
        new_clus[clu_index][i] += sub_mat[i]

print('sum_dis:',sum_dis)
# print new centers
for i in range(CLU_NUM):
    if new_clu_sums[i] < 10000:
        continue
#     print(new_clu_sums[i],end=' ')
    for j in range(PATTEN_WID):
        print(new_clus[i][j]/new_clu_sums[i],end=' ')
    print()
