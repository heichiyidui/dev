#!/usr/bin/env python3
import numpy as np

#######################################
# load the Pi array and the Q matrix from the last round of estimation
Pi = np.genfromtxt('vtml/Pi')
F = np.diag(Pi)
Q = np.genfromtxt('vtml/Q')

ev = np.linalg.eig(Q)
Lambda = ev[0]
S = ev[1]
inv_S = np.linalg.inv(S)

#######################################
# read the alignments and find the alignment distances
AA_TO_INT = {'A': 0 ,'R': 1 ,'N': 2 ,'D': 3 ,'C': 4 ,\
             'Q': 5 ,'E': 6 ,'G': 7 ,'H': 8 ,'I': 9 ,\
             'L': 10,'K': 11,'M': 12,'F': 13,'P': 14,\
             'S': 15,'T': 16,'W': 17,'Y': 18,'V': 19,};

import scipy.optimize

def loglikeN(t): # log likelihood of time t
# in R:    -sum(n_mat * log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)) )
    diag_m = np.diag(np.exp(Lambda*t))
    a_mat = np.dot(S,diag_m)
    b_mat = np.dot(a_mat,inv_S)
    c_mat = np.dot(F,b_mat)
    d_mat = np.multiply(aln_mat,np.log(c_mat))
    return -np.sum(np.multiply(aln_mat,np.log(c_mat)))

import sys
dom_ls = open(sys.argv[1]).read().split()
# dom_ls = open('index/dom.ls').read().split()

for dom_id in dom_ls:
    print('>'+dom_id)
    ifile=open('../cath/bl_out/'+dom_id)
    ifile.readline()
    dom_seq = ifile.readline()[:-1]
    for line in ifile:
        aln_seq = ifile.readline()[:-1]
        aln_mat = np.outer(Pi,Pi)       # Laplace smoothing?
        for i in range(len(aln_seq)):
            if aln_seq[i] == '-':
                continue
            c_i = AA_TO_INT[dom_seq[i]]
            c_j = AA_TO_INT[aln_seq[i]]
            aln_mat[c_i][c_j] += 1
            aln_mat[c_j][c_i] += 1

        aln_mat /= np.sum(aln_mat)
        print(scipy.optimize.minimize(loglikeN,0.1).x[0])

    ifile.close()
