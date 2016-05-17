#!/usr/bin/env python3

import numpy as np

ALPHA = 0.02
n_sum_mat = 84

#######################################
# load the sum matrics
Pi = np.genfromtxt('cao/Pi_5')
Pi.shape = (400)

sum_aln_mats=[]
for i in range(100):
    in_file_name = 'sum_mat/{:03d}'.format(i)
    sum_aln_mats.append( np.genfromtxt(in_file_name) )

p_mats = []
for i in range(100):
    p_mat = sum_aln_mats[i]
    for m in range(400):
        p_mat[m] /= sum(p_mat[m])
    p_mats.append(p_mat)

#######################################
#  get the resolvent

cao_dis_strs = open('cao.dis').read().split()
cao_dis = [float(x) for x in cao_dis_strs]
cao_dis_ranges = []
cao_dis_ranges.append((cao_dis[0]+cao_dis[1])/2)
for i in range(1,n_sum_mat):
    t_range = (cao_dis[i+1] - cao_dis[i-1]) / 2
    cao_dis_ranges.append(t_range)

R = np.zeros((400,400),dtype=float)
for i in range(n_sum_mat):
    R += np.exp(-ALPHA*cao_dis[i]) * p_mats[i] * cao_dis_ranges[i]

inv_R = np.linalg.inv(R)
Q = ALPHA * np.eye(400) - inv_R

#######################################
# normalize Q

# Q should be semi-symmetric
for i in range(400):
    Q[i] = Q[i] * Pi[i]

Q = (Q + Q.T) * 0.5

for i in range(400):
    Q[i] = Q[i] / Pi[i]

# off-diag elements should all be non-negtive
for i in range(400):
    for j in range(400):
        if i==j :continue
        if Q[i,j] < 0: Q[i,j]=0

# sum(Q.T) should be zero
t_mat = Q - np.diag(np.diag(Q))
Q = t_mat - np.diag(sum(t_mat.T))

# I-Q should give PAM1

Q /= np.sum(t_mat) / 4 # 0.01 * 400

#######################################
# write the Q matrix
np.savetxt('cao/Q',Q,'%08e')

