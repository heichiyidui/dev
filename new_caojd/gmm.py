#!/usr/bin/python2.7

import itertools
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn import mixture

X=[]
ifile=open('t.in')
for line in ifile:
    cols=line.split()
    x=[float(cols[0]),float(cols[1]),float(cols[2])]
    X.append(x)

ifile.close()
X=np.array(X)

# lower_bounds=[]
lowest_bic = np.infty
bic = []
n_components_range = range(2, 18)
for n_components in n_components_range:
    # Fit a mixture of Gaussians with EM
    gmm = mixture.GMM(n_components=n_components, \
                        covariance_type='diag',n_iter=500)
    gmm.fit(X)
    bic.append(gmm.bic(X))
    # z=best_gmm.predict_proba(X)
    # lower_bounds.append(best_gmm.lower_bound(X,z))
    if bic[-1] < lowest_bic:
        lowest_bic = bic[-1]
        best_gmm = gmm

bic = np.array(bic)
best_gmm.weights_
best_gmm.means_
# best_gmm.precs_
best_gmm.converged_

# lower_bounds

np.savetxt('t.out',best_gmm.predict_proba(X))
gmm=best_gmm
np.savetxt('gmm.means',gmm.means_)
np.savetxt('gmm.weights',gmm.weights_)
np.savetxt('gmm.covars',gmm.covars_)