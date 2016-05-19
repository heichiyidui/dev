#!/usr/bin/env python3
import numpy as np

#######################################
# load the Pi array and the Q matrix
Pi = np.genfromtxt('cao/Pi_5')
Pi.shape = (400)

F = np.diag(Pi)
Q = np.genfromtxt('cao/Q')

ev = np.linalg.eig(Q)
Lambda = ev[0]
S = ev[1]
inv_S = np.linalg.inv(S)

import sys
t = int(sys.argv[1])

diag_m = np.diag(np.exp(Lambda*t))
a_mat = np.dot(S,diag_m)
b_mat = np.dot(a_mat,inv_S)
c_mat = np.dot(F,b_mat)

base_mat = np.outer(Pi,Pi)

score_mat = np.log(c_mat/base_mat)
score_mat = (score_mat + score_mat.T)/2

np.savetxt('t.out',score_mat,fmt='%7.4f')
