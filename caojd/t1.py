#!/usr/bin/env python3

MAX_DIS = 241
ALPHA = 0.01

import numpy as np

#######################################
# 3. normalize the sum matrices, write them to files
sum_aln_mats=[]
for i in range(100):
    ifile_name = 'sum_mat/{:03d}'.format(i)
    sum_aln_mats.append(np.genfromtxt(ifile_name))

#######################################
# 4. get the P(t) matrices

p_mats = []
for i in range(100):
    p_mat = sum_aln_mats[i]
    for m in range(20):
        p_mat[m] /= sum(p_mat[m])
    p_mats.append(p_mat)

#######################################
# 5. get the resolvent

dis_step = MAX_DIS/100
half_step = dis_step/2

dis = []
for i in range(100):
    dis.append(MAX_DIS*i/100 + half_step)

R = np.zeros((20,20),dtype=float)
for i in range(100):
    R += np.exp(-ALPHA*dis[i]) * p_mats[i] * dis_step

inv_R = np.linalg.inv(R)

Q = ALPHA * np.eye(20) - inv_R

#######################################
# 6. normalize Q

# sum(Q) should be zero

Q -= np.ones((20,20),dtype=float) * np.sum(Q) /400

# I-Q should give PAM1

t_mat = Q - np.diag(np.diag(Q))

Q /= np.sum(t_mat) / 0.2

#######################################
# 7. write the Q matrix
np.savetxt('vtml/Q',Q,'%08e')



