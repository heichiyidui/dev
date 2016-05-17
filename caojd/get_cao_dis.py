#!/usr/bin/env python3

import numpy as np

#######################################
# load Pi and the sum matrics
Pi = np.genfromtxt('cao/Pi_5')
Pi.shape = (400)
F = np.diag(Pi)
Q = np.genfromtxt('cao/Q')

sum_aln_mats=[]
for i in range(100):
    in_file_name = 'sum_mat/{:03d}'.format(i)
    sum_aln_mats.append( np.genfromtxt(in_file_name) )

#######################################

ev = np.linalg.eig(Q)
Lambda = ev[0]
S = ev[1]
inv_S = np.linalg.inv(S)

import scipy.optimize

def loglikeN(t): # log likelihood of time t
# in R:    -sum(n_mat * log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)) )
    diag_m = np.diag(np.exp(Lambda*t))
    a_mat = np.dot(S,diag_m)
    b_mat = np.dot(a_mat,inv_S)
    c_mat = np.dot(F,b_mat)
    d_mat = np.multiply(aln_mat,np.log(c_mat))
    return -np.sum(d_mat)

for aln_mat in sum_aln_mats:
    print(scipy.optimize.minimize(loglikeN,0.1).x[0])
