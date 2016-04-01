################################################################################
#                               YASPIN                                         #
#                                                                              #
#    Yet Another Secondary structure Predictor usIng Neual network             #
################################################################################

################################################################################
#                                                                              #
# The method was published in Bioinformatics. 21(2):152-9. (Epub 2004 Sep 17). #
#                                                                              #
# It uses a hidden neural network model for protein secondary structure        #
# prediction.                                                                  #
#                                                                              #
# Three-states secondary structures are re-classified to 7 states.             #
# Helix begins and ends, strand begins and ends are seperated from the helix   #
# and stand classes.                                                           #
#                                                                              #
# A simple back-propagation supervised feed-forward neural network was trained #
# for the prediction of emission probabilities of the 7 states using           #
# PSSMs (Position Specific Scoring Matrices) from PsiBlast search.             #
#                                                                              #
# A hidden Markov model is then applied on the emission probabilities to       #
# assign the likely transimissions between states.                             #
#                                                                              #
################################################################################

################################################################################
# I think now it is about the time to upgrade the method.                      #
################################################################################

################################################################################
#                                                                              #
#  Some ideas to be tested...                                                  #
#                                                                              #
#  1. more than 7 states? beta-bridge to be seperated from strand?             #
#     3/10 helix to be seperated from helix? pi helix, hydrogen bonded turn    #
#     and bend to be used? helix begins and ends are to be included into the   #
#     helix state?                                                             #
#                                                                              #
#  2. The consensus alignment to be used to find and fill the unknown gaps in  #
#     query sequences. (implemented in the CATH scripts)                       #
#                                                                              #
#  3. A second order hidden Markov model with both secondary structure and     #
#     burial states?                                                           #
#                                                                              #
#  4. The snow ball approach for neural network training?                      #
#                                                                              #
#  5. Include prediction from matching to known structures?                    #
#                                                                              #
################################################################################

################################################################################
# 1. select domain set                                                         #
################################################################################

# from CATH v3.5.0, we obtained 41213 QCed domains from the CATH database.
# 11221 S35 families.

# from CATH v4.0.0, we obtained 16556 QCed S35 family representative domains. 

# 7/8 states secondary structure definitions by DSSP and STRIDE 
# are somewhat different (20 +- 7 % difference)
# 3 states definitions are highly similar (5 +- 3% difference)
# removed 104 domains with > 15% differences 

mkdir index 
index/s.ls 
index/test.ls 
index/train.ls 

# 16452 S35 family representative domains
# 481 (2.8%) removed from the original 16933 
# considering 275 were removed as transmembrane, no bad at all. 

# 1/3 (5484)  as the test set 
# 2/3 (10968) as the training set 

################################################################################
# 2. Statistics of the residue states                                          #
################################################################################

# DSSP definitions are well accepted. 
# STRIDE defintions seem to be better.
# (not on the few domains with very different DSSP and STRIDE defs.)
# I'm using STRIDE ACC, so I'm also using STRIDE SS. 
# It's also easier because STRIDE assigned all residues, DSSP missed a few.

# 2442984 residues in 16452 s35 domains 
# 
# H 944554    0.3866
# E 564943    0.2313
# C 933487    0.3821
#
#
# H 851809    0.3487
# G  92745    0.0380
# E 540229    0.2211
# B  24714    0.0101
# C 445826    0.1825
# T 487336    0.1995
# I    325    0.0001
# S      0    0.0

get_rsa.py > t.in 
listdis t.in -bn 781 > t.dat 
# binnum:    781
tail -n +2 t.dat > t2.dat 

# RAC
# min:       0
# max:       2.3837
# mean:      0.339529
# median:    0.302206 
# sdv:       0.28418

# 7.5% are 0 

# using R package mixdist
# treat the observed RAC as a mixture of two distributions

R

library(mixdist)
rac_data=as.matrix(read.table('t2.dat',header=FALSE))
fit=mix(rac_data,mixparam(c(0.1,0.6),.5), "weibull")
# two components mixture wouldn't fit the data well 
# fit$arameters:
#       pi     mu  sigma
# 1 0.4042 0.1735 0.2085
# 2 0.5958 0.5151 0.2230

fit=mix(rac_data,mixparam(c(0.1,0.2,0.6),.5),"weibull")
# three components mixture looks better
# but this is only a local optimal solution

fit=mix(rac_data,mixparam(c(0.05,0.3,0.6),.5),"weibull")
# this will get the global optimal fit
# fit$parameters
#          pi        mu     sigma
# 1 0.3378852 0.1440214 0.2050568
# 2 0.1152352 0.2320332 0.1267487
# 3 0.5468796 0.5385876 0.2176676

weibullpar(0.1440214,0.2050568)
#      shape     scale loc
#   0.716621 0.1162887   0

weibullpar(0.2320332,0.1267487)
#      shape     scale loc
#   1.904627 0.2615063   0

weibullpar(0.5385876,0.2176676)
#      shape     scale loc
#   2.664756 0.6059084   0

#    pi1 = 0.3378852 ; k1 = 0.716621 ; l1 = 0.1162887 
#    pi2 = 0.1152352 ; k2 = 1.904627 ; l2 = 0.2615063 
#    pi3 = 0.5468796 ; k3 = 2.664756 ; l3 = 0.6059084 

#    
#    p1 = pi1 * (k1/l1) * ((x/l1) ** (k1-1)) *exp(- (x/l1)**k1 )
#    p2 = pi2 * (k2/l2) * ((x/l2) ** (k2-1)) *exp(- (x/l2)**k2 )
#    p3 = pi3 * (k3/l3) * ((x/l3) ** (k3-1)) *exp(- (x/l3)**k3 )

#    p=  p1 + p2 + p3 

# use p1 + p2 as buried
# cut at 27.1% 
# 45.3% as buried, 54.7% as exposed. 

# remember to set p1(0)=1, p2(0)=0 and p3(0)=0
# and p1(x>0.75)=0, p2(x>0.75)=0, p3(x>0.75)=1

get_exp.py > index/cath_s35.exp 

# the plot of mixtures
res/exp_fit.png 

################################################################################
# 3. Neural Network implementation                                             #
################################################################################

yaspin_1.0_src/
# My old C++ implementation of network, sigmoid transition for the hidden layer,
# softmax for the output layer, supervised training. 
# Found in the old YASPIN package. 
# Its performance on the test set of 812 domains, 109031 residues (CATH v3.5),
# is pretty good...

snn.py 
# The Python implementation is way shorter.
# A simple feedforward back-propagation neural network in Python 
# It has one hidden layer, with the tanh activity function.
# The output layer uses the softmax activity function.
# copied a lot from Jose Antonio Martin 
# http://www.dacya.ucm.es/Jam/download.htm

# 1630993 residues in the training set 
#  residues in the test set 

cv_expos_train.py
cv_expos_train.sh
for hu in 10 20 30 40 60 80 100; do 
    for rs in 0 1 2 3 4 5 6 7 8 9; do 
        qsub cv_expos_train.sh $hu $rs 
    done
done 

# test the networks
expos_nn_test.py

# best cv_error 0.4982 from a network of 30 hidden units, continous targets
# best accuracy 75.96% from a network of 20 hidden units, binary targets

# 69~76% are correctly predicted.
# The cross validation (corss entropy) error is anti-correlated to the 
# prediction accuracy. The correlation coefficient is -0.6732289.
# The t-value of the correlatioin coefficient is -9.890182. It is significant.

# No obvious difference between using binary or continuous training targets.
# No obvious difference between different hidden layer sizes.
# Smaller networks are better

# big differences between the training results from the networks with the same 
# layout but different random seeds.

# Try momentum factor of 0.5 instead of 0.8 
# with learning rate from 0.01 to 0.001
# will NOT converge in 1000 iteartions

# try [0.8, 0.0001]
for hu in 20 30 40 ; do 
    for rs in 0 1 2 3 4 5 6 7 8 9; do 
        qsub cv_expos_train.sh $hu $rs 
    done
done 

# will converge in 100~800 iterations.
# last 1~7 days
# cross validation error 0.552~0.559
# no much difference given different hidden layer sizes. 

# try [0.7, 0.001]
# will converge soon after the first 30 'burn in' iterations
# cross validation error 0.552~0.619 

# try [0.8, 0.0005]
# will converge in 33 ~ 75 iterations
# cross validation error 0.550~0.594
# more iterations means smaller errors

# try [0.9, 0.0002]
# converge in 33 ~ 45 iterations 
# error 0.556 ~ 0.562

# switch to 10 iterations for 'burning in'
# retry [0.9, 0.0002]
# converge in 15 ~ 51 iterations
# error 0.556 ~ 0.584

# will stay with [0.9, 0.0002], 20 hidden units

################################################################################
# 4. Secondary structure states definitions                                    #
################################################################################

#######################################
# 4.1 8/7 states STRIDE definitions 

#  24848 B 1 
# 543272 E 0.04574
#  94313 G 0.26346
# 853568 H 0.02911
#    340 I
# 492100 T 0.05049
# 448374 C 0.05542
# 239253 X

# 'I' is too few. 
# Can we tell 'B' from 'E', 'G' from 'H' or 'T' from 'C'?
# use about 24848 residues for each state, 6-fold cross-validation 

test_8_state.py -s1 C -s2 H -rs 0 > C_H_0.out 
test_8_state.py -s1 C -s2 H -rs 1 > C_H_1.out 
# ...
test_8_state.py -s1 C -s2 E -rs 0 > C_E_0.out 
# ...
test_8_state.py -s1 T -s2 C -rs 9 > T_C_9.out 

# 'I' is too few, 'S' is missing in STRIDE definitions
# We end up dealing with 6 states.
# using cross-validation errors to test the predictablity between states. 
# 1.000 for the same state

#       B      E      G      H      C      T
# B     1.0000 0.4684 0.4776 0.4776 0.5937 0.5632
# E     0.4684 1.0000 0.3440 0.2547 0.4441 0.3827
# G     0.4776 0.3440 1.0000 0.4725 0.5137 0.5676
# H     0.4776 0.2547 0.4725 1.0000 0.3576 0.3420
# C     0.5937 0.4441 0.5137 0.3576 1.0000 0.6224
# T     0.5632 0.3827 0.5676 0.3420 0.6224 1.0000

# 'E' is different from most states but 'B'.
# 'B' is too few. 
# Merge 'E' and 'B'.

# left with 5 states, 'E', 'G', 'H', 'T', 'C'

# 448374 C 0.21034
# 568120 E 0.16601
#  94313 G 1.00000
# 853568 H 0.11049
# 492100 T 0.19165 

test_5_state.py -s1 C -s2 H -rs 0 > C_H_0.out
# ...

#   E       G       H       C       T 
# E 1.0000  0.3260  0.2286  0.4672  0.4036
# G 0.3260  1.0000  0.4617  0.4953  0.5769
# H 0.2286  0.4617  1.0000  0.3294  0.3342
# C 0.4672  0.4953  0.3294  1.0000  0.6356
# T 0.4036  0.5769  0.3342  0.6356  1.0000

# still can't decide whether to seperate 'G' and 'H' or not 

# Going DSSP from now on. At least it doesn't have that much 
# (64 in stead of 1138) single residue 'G'

#######################################
# 4.2 8/6 states DSSP definitions 

#  25813 B
# 527956 C
# 518361 E
#  88085 G
# 814885 H
#    457 I
# 207374 S
# 273774 T

# 'I' too few, treated as 'C'
# 'B' too few, treated as 'E'

# 528413 C 0.16670
# 544174 E 0.16187
#  88085 G 1.00000
# 814885 H 0.10810
# 207374 S 0.42476
# 273774 T 0.32174

# try to test 'G' vs 'H' and 'S' vs 'C' and 'T' vs 'C'

test_6_state.py -s1 C -s2 H -rs 0 > C_H_0.out
# ...

# 150 one hour jobs on the cluster 

#     E      G      H      T      S      C 
# E   1.000  0.326  0.227  0.319  0.398  0.477
# G   0.326  1.000  0.412  0.576  0.534  0.521
# H   0.227  0.412  1.000  0.369  0.302  0.339
# T   0.319  0.576  0.369  1.000  0.592  0.527
# S   0.398  0.534  0.302  0.592  1.000  0.643
# C   0.477  0.521  0.339  0.527  0.643  1.000

# 'T' and 'S' can not be seperated from 'C'.
# 'G' is still trouble. It can be seperated from 'H'.
# However, the network can not efficiently identify it from 'T', 'S' and 'C'. 
# The established secondary structure classes 'E' and 'H' are not like that. 

# Just treat 'G' as 'H'.

# 'E','B'         -> 'E'
# 'H','G'         -> 'H'
# 'C','S','T','I' -> 'C'

# 3 states
# 1009561 C  41.09 %
#  544174 E  22.15 %
#  902970 H  36.76 %

# to get the 7 states definitions:
# 'HEE','CEE' -> 'D'
# 'EEH','EEC' -> 'F'
# 'CHH','EHH' -> 'G'
# 'HHC','HHE' -> 'I'

# 1009561 C 0.09439
#   97818 D 0.97419
#  348541 E 0.27341
#   97815 F 0.97422
#   95353 G 0.99937
#  712324 H 0.13378
#   95293 I 1.00000

#######################################
# 4.3 to test in the 7 states definitions if we can merge 'D', 'E' and 'F'

test_7_state.py -s1 C -s2 H -rs 0 > C_H_0.out
# ...

#    C      D      E      F      G      H      I 
# C  0.5000 0.4004 0.3872 0.4281 0.4085 0.3357 0.4469
#
# D  0.4004 0.5000 0.4779 0.2675 0.2143 0.1994 0.2383
# E  0.3872 0.4779 0.5000 0.4979 0.2109 0.2017 0.2352
# F  0.4281 0.2675 0.4979 0.5000 0.2210 0.1813 0.2312
#
# G  0.4085 0.2143 0.2109 0.2210 0.5000 0.3576 0.2668
# H  0.3357 0.1994 0.2017 0.1813 0.3576 0.5000 0.4230
# I  0.4469 0.2383 0.2352 0.2312 0.2668 0.4230 0.5000

# very tempting to merge 'D', 'E' and 'F'

################################################################################
# 5. Secondary structure prediction                                            #
################################################################################

#######################################
# 5.1 basic 7-states prediction

# to get the 7-state definitions:
dssp_to_seven.py 

# train it with cross-validation 
train_7_state.py -hu 30 -rs 2

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m e
#$ -N y7_train
#$ -cwd
#$ -q long.q,bignode.q
#$-l h_vmem=3G

hu=$1
rs=$2

./train_7_state.py -hu $hu -rs $rs  > y7_tra_${hu}_${rs}.out
#--------------------------------------

for hu in 10 20 30 40 60 80 100 120 150 ; do 
    for rs in 0 1 2 3 4 5 6 7 8 9; do 
        qsub t.sh $hu $rs 
    done
done 

# the best network is y_7_80_1.net 
res/y_7_80_1.net
# cross-validation error 0.85398
# test error 0.89075

# converged after 84 iterations. 
# 60-80 hidden units should be optimal. 
# training time is about 6 hours to 3 days.
# most should be fininshed in about 1 day. 

#######################################
# 5.2 balanced data set

# now, we know the training set is biased. 
# 'G' and 'I', maybe 'D' and 'F' states are about 10 time less. 
# try randomly select subsets for 'C', 'E' and 'H' so the states presented even.
# under sampling for a balanced training set 
bal_tra_7_state.py

res/y_7_bal_80_6.net
# cross-validation error 1.10049 on the balanced set.
# test error 1.26700 on the unbalanced test set.
# 80 hidden units

test_y7_net.py res/y_7_80_1.net

test_y7_net.py res/y_7_bal_80_6.net

# the mean CE error is 1.51103 if we use the background distribution as output
#                      1.94591 if use [1/7] * 7 as output

# unbalanced training
# CE error 0.89075
#       real_states  avg_out  avg_real_out
# 332512 C  0.40952  0.37865  0.62513
#  32229 D  0.03969  0.04379  0.25092
# 114946 E  0.14157  0.15615  0.53103
#  32226 F  0.03969  0.04674  0.23346
#  31709 G  0.03905  0.03321  0.31442
# 236656 H  0.29146  0.30250  0.73793
#  31680 I  0.03902  0.03895  0.18721
# sum avg_real_out is 2.88010
# 
# balanced training
# CE error 1.26700
#       real_states  avg_out  avg_real_out
# 332512 C  0.40952  0.19971  0.33690
#  32229 D  0.03969  0.10489  0.48042
# 114946 E  0.14157  0.12849  0.42078
#  32226 F  0.03969  0.11001  0.43099
#  31709 G  0.03905  0.12149  0.55009
# 236656 H  0.29146  0.18322  0.47919
#  31680 I  0.03902  0.15220  0.48889
# sum avg_real_out is 3.18726
# 

# balanced training followed by unbalanced training
# CE training error 0.86929
# CE error 0.93293
#       real_states  avg_out  avg_real_out
# 332512 C  0.40952  0.33774  0.52617
#  32229 D  0.03969  0.03572  0.17877
# 114946 E  0.14157  0.22856  0.62691
#  32226 F  0.03969  0.03766  0.15906
#  31709 G  0.03905  0.03971  0.32894
# 236656 H  0.29146  0.28130  0.68292
#  31680 I  0.03902  0.03931  0.17230
# sum avg_real_out is 2.67507

# balanced (via under-sampling) training doesn't work for multi-class targets?
# network ensemble? thresholding?

#######################################
# 5.3 network ensemble 

test_y7_net_ensemble.py 
# imbalance network ensemble:
# CE error 0.83340
#       real_states  avg_out  avg_real_out
# 332512 C  0.40952  0.39031  0.62617
#  32229 D  0.03969  0.04182  0.23614
# 114946 E  0.14157  0.14393  0.49710
#  32226 F  0.03969  0.04154  0.19924
#  31709 G  0.03905  0.03884  0.31432
# 236656 H  0.29146  0.30539  0.72206
#  31680 I  0.03902  0.03817  0.16999
# 
# some small improvements over the single network
# CE is much better
# the main difference is the reduction of false-negtives 
#
# balanced network ensemble:
# CE error 1.18714
#       real_states  avg_out  avg_real_out
# 332512 C  0.40952  0.17998  0.30806
#  32229 D  0.03969  0.11568  0.48986
# 114946 E  0.14157  0.11912  0.39280
#  32226 F  0.03969  0.11963  0.46113
#  31709 G  0.03905  0.12247  0.55053
# 236656 H  0.29146  0.20282  0.51774
#  31680 I  0.03902  0.14031  0.47915
# again, improve the result by reducing the false-negtives
 res/net_ensemble.png

# the best y7_net here is res/y7_nets/y_7_80_47.net
# with cross-validation CE error 0.84825
# test error 0.87974

#######################################
# 5.4 prediction confidence vs error 

# Define confidence as the Kullback–Leibler divergence between the background 
# distribution of the states and the prediction. 

# from numpy import * 
# back_ground=array([0.40952,0.03969,0.14157,0.03969,0.03905,0.29146,0.03902])
# confidence = dot(log(out/back_ground),out)

# Linear regresion between CE error and confidence.
# Overall, the the correlation coefficient is negtive. 
# When the network is more confident, the prediction error is small.

# When the real state is 'C', the correlation coefficient is positive. 
# With other states, it is always negtive. 

# It is possible to sum networks outputs with their confidences as weights.

####################################
# 5.6 conclusions

# Some interesting observations:
# 1. network ensembles are better, especially on reducing false-negtives.
# 2. with imbalance data, H, C, E and G can be predicted. 
# 3. with balanced (under-sampled) data, prediction of C and E suffers. 
#    G and H can still be predicted. 


# Conclusion:
# 1. must use network ensemble
# 2. do not under-sample
# 3. use H, C, E, G and H-begin states

################################################################################
# 6. HMM model                                                                 #
################################################################################

########################################
# 6.1 the old 7 states HMM model 
# 
# the prior distribution of 7 states:
#
# 1003641  C  0  0.41084
#   97256  D  1  0.03981
#  346650  E  2  0.14190
#   97259  F  3  0.03981
#   94822  G  4  0.03882
#  708482  H  5  0.29002
#   94765  I  6  0.03879
#
# The forward and backward transition tables:
# obviously they are the transposed forms to each other
#        C       D       E       F       G       H       I
# C 767818   94672   21633       0   88594      77       0
# D      0       0   83695   13561       0       0       0
# E  20942       0  239788   83698    1729       1       0
# F  92757       0       0       0    4499       3       0
# G      0       0       0       0       0   94822       0
# H     86       5       1       0       0  613021   94765
# I  91154    2579    1032       0       0       0       0
#
 
#        C       D       E       F       G       H       I
# C 767818       0   20942   92757       0      86   91154 
# D  94672       0       0       0       0       5    2579
# E  21633   83695  239788       0       0       1    1032
# F      0   13561   83698       0       0       0       0
# G  88594       0    1729    4499       0       0       0
# H     77       0       1       3   94822  613021       0
# I      0       0       0       0       0   94765       0
#
# The forward and backward transition probabilities tables:
# back_ground=[0.41084, 0.03981, 0.14190, 0.03981, 0.03882, 0.29002, 0.03879]

# 0.78929107,0.09731962,0.02223813,0.00000004,0.09107164,0.00007945,0.00000004
# 0.00000422,0.00000041,0.86055648,0.13943510,0.00000040,0.00000298,0.00000040
# 0.06049940,0.00000012,0.69271099,0.24179074,0.00499493,0.00000373,0.00000011
# 0.95370564,0.00000041,0.00000146,0.00000041,0.04625785,0.00003383,0.00000040
# 0.00000433,0.00000042,0.00000150,0.00000042,0.00000041,0.99999251,0.00000041
# 0.00012207,0.00000712,0.00000161,0.00000006,0.00000005,0.86599728,0.13387180
# 0.96188940,0.02721482,0.01089148,0.00000042,0.00000041,0.00000306,0.00000041

# 0.78932110,0.00000004,0.02152862,0.09535469,0.00000004,0.00008871,0.09370680
# 0.97342516,0.00000041,0.00000146,0.00000041,0.00000040,0.00005439,0.02651777
# 0.06249721,0.24178836,0.69272900,0.00000012,0.00000011,0.00000373,0.00298148
# 0.00000422,0.13943080,0.86056078,0.00000041,0.00000040,0.00000298,0.00000040
# 0.93431352,0.00000042,0.01823547,0.04744671,0.00000041,0.00000306,0.00000041
# 0.00010935,0.00000006,0.00000161,0.00000429,0.13394362,0.86594101,0.00000005
# 0.00000434,0.00000042,0.00000150,0.00000042,0.00000041,0.99999251,0.00000041


########################################
# 6.2 the implementation of the forward-backward algorithm
#

# 6.2.1 on the old y7 structure
test_yaspin_7.py 
# 5484 domains, 891095 residues, 811958 residues with structure definitions
# 77.78 % correct 3-states prediction
# 70.96 % correct 7-states prediction 
#     real     pred 
# 0 332512   354004
# 1  32229    22284
# 2 114946   135735
# 3  32226    20940
# 4  31709    18539
# 5 236656   245086
# 6  31680    15370

# marginal likelihood ranges from 0.226778 to 1
# in 10 bins:
#                          # residues    Q3
# 0.0000000   0.6735415    81196         54.43 %
# 0.6735415   0.8498510    81195         60.17 %
# 0.8498510   0.9477252    81197         65.32 %
# 0.9477252   0.9857625    81196         71.32 %
# 0.9857625   0.9966643    81195         77.07 %
# 0.9966643   0.9992816    81199         82.24 %
# 0.9992816   0.9998613    81201         86.59 %
# 0.9998613   0.9999798    81119         90.15 %
# 0.9999798   0.9999990    82137         93.68 %
# 0.9999990   1.0000000    80323         96.95 %

# 6.2.2 network ensemble
test_yaspin_7_esb.py
# 78.42 % correct 3-states prediction
# 71.50 % correct 7-states prediction 
#     real     pred 
# 0 332512   340557
# 1  32229    19756
# 2 114946   141987
# 3  32226    15906
# 4  31709    18652
# 5 236656   262470
# 6  31680    12630

# max marginal likelihood ranges from 0.23386 to 1 
# in 10 bins 
#                          # of residue  Q3
# 0.0000000   0.6387860    81195         54.30 %
# 0.6387860   0.8110115    81196         61.06 %
# 0.8110115   0.9248632    81197         66.20 %
# 0.9248632   0.9761472    81194         71.96 %
# 0.9761472   0.9932850    81197         77.93 %
# 0.9932850   0.9982654    81195         83.25 %
# 0.9982654   0.9996067    81195         87.44 %
# 0.9996067   0.9999347    81176         90.93 %
# 0.9999347   0.9999979    81401         94.21 %
# 0.9999979   1.0000000    81012         96.92 %

# Don't know why there's underprediction of the minor states. 

################################################################################
# 7. new yaspin with 10 states                                                 #
################################################################################

#######################################
# 7.1 network training 

# with 5 secondary structure states: C, E , Hb, H , G
# with 2 exposure states: buried and exposed. 

# 10 states
#          C  E Hb  H  G
# exposed  0  1  2  3  4
# buried   5  6  7  8  9

# in 16452 S domains
# 1003641 C
#  541165 E
#   87583 G
#  739489 H
#   70997 Hb

#  661718 0   341923 5
#  174656 1   366509 6
#   40690 2    30307 7
#  366806 3   372683 8
#   55377 4    32206 9

# 'C' and 'G' are more often exposed.
# 'Hb' are a bit more exposed. 
# 'E' are more often buried. 

train_10_state.py -hu 80 -rs 9

# 100 networks trained. 
# CE cross-validation error ranges between 1.18523 and 1.39030
# best network y_10_80_02.net

test_y10_net.py
# test CE error  1.23502
#    real     avg_out  avg_real_out
# C  0.27088  0.25374  0.46685
# E  0.07150  0.10135  0.33584
# Hb 0.01666  0.01962  0.31012
# H  0.15015  0.14468  0.49913
# G  0.02267  0.01435  0.05774
# C  0.13997  0.10585  0.22153
# E  0.15003  0.17718  0.55666
# Hb 0.01241  0.01103  0.16291
# H  0.15256  0.16352  0.58189
# G  0.01318  0.00869  0.03851

# looks like we have serious under-prediction of 'G' states. 

test_y10_net_ensemble.py

# test CE error  1.14716
#    real     avg_out  avg_real_out
# C  0.27088  0.28056  0.52316
# E  0.07150  0.07575  0.27310
# Hb 0.01666  0.01650  0.26042
# H  0.15015  0.15618  0.51019
# G  0.02267  0.02183  0.08017
# C  0.13997  0.12793  0.27434
# E  0.15003  0.14524  0.50843
# Hb 0.01241  0.01274  0.19010
# H  0.15256  0.15068  0.53471
# G  0.01318  0.01260  0.05138

# under-prediction of 'G' and 'Hb' was improved...
#######################################
# 7.2 transition table 

get_trans_table.py

# The first order transition tables are easy to get. 

# forward transition table
#         0      1      2      3      4      5      6      7      8      9  
# 0 [367955, 28974, 28120,    24,  7745,150169, 33090, 15290,    11,  4158,],
# 1 [ 22013, 55282,   218,     1,   314, 18954, 75884,  1017,     0,   717,],
# 2 [     0,     0,     0, 32661,     0,     0,     0,     0,  8029,     0,],
# 3 [ 30606,    28,     0,166285,   230, 14889,   135,     0,154021,   309,],
# 4 [  5236,    88,     0,   234, 27941,  7312,   842,     0,   869, 12761,],
# 5 [157114, 18388, 11908,    10,  5173, 92580, 35853, 12159,     4,  4069,],
# 6 [ 33468, 71312,   444,     1,   614, 39264,218264,  1841,     0,  1065,],
# 7 [     0,     0,     0, 18060,     0,     0,     0,     0, 12247,     0,],
# 8 [ 15039,    92,     0,148535,   973, 10433,   525,     0,196156,   750,],
# 9 [  3987,   242,     0,   691, 12311,  3738,  1665,     0,  1206,  8339,],

# backward transition table 
#         0      1      2      3      4      5      6      7      8      9  
# 0 [367955, 22013,     0, 30606,  5236,157114, 33468,     0, 15039,  3987,],
# 1 [ 28974, 55282,     0,    28,    88, 18388, 71312,     0,    92,   242,],
# 2 [ 28120,   218,     0,     0,     0, 11908,   444,     0,     0,     0,],
# 3 [    24,     1, 32661,166285,   234,    10,     1, 18060,148535,   691,],
# 4 [  7745,   314,     0,   230, 27941,  5173,   614,     0,   973, 12311,],
# 5 [150169, 18954,     0, 14889,  7312, 92580, 39264,     0, 10433,  3738,],
# 6 [ 33090, 75884,     0,   135,   842, 35853,218264,     0,   525,  1665,],
# 7 [ 15290,  1017,     0,     0,     0, 12159,  1841,     0,     0,     0,],
# 8 [    11,     0,  8029,154021,   869,     4,     0, 12247,196156,  1206,],
# 9 [  4158,   717,     0,   309, 12761,  4069,  1065,     0,   750,  8339,],

# The third order transition table got many empty columns. 
# 640 in the forward table, 628 in the backward table.
# Some combinations of 3 states were never observed. 

# back_ground=[0.270877, 0.071496, 0.016657, 0.150153, 0.022669,
#              0.139967, 0.150032, 0.012406, 0.152559, 0.013184]

# worrying Hb and G are too few...

# use the back_ground table for these columns

#######################################
# 7.3 10-states first-order HMM implementation

test_yaspin_10.py
# Q3  74.72 %
# Q10 56.91 %

#    real    predicted  true-positive
# 0  220065  237988     142469    64.740 %
# 1   58052   71612      26664    45.931 %
# 2   13708   12258       5930    43.259 %
# 3  123002  116494      76171    61.927 %
# 4   18364    4467       1100    05.990 %
# 5  112447   39805      20596    18.316 %
# 6  121349  178118      94040    77.495 %
# 7   10096    4820       2251    22.296 %
# 8  124313  144164      92349    74.287 %
# 9   10562    2232        489    04.630 %

# serious under-prediction of 'G' states
# 'Hb' exposed is fine, but not 'Hb' buried. 
# 'C'  exposed is under-predicted.

test_yaspin_10_esb.py

# Q3  77.94 %
# Q10 59.79 %

#    real    predicted  true-positive
# 0  220065  270174     160032    72.720 %
# 1   58052   42746      20861    35.935 %
# 2   13708    9853       5461    39.838 %
# 3  123002  123771      80570    65.503 %
# 4   18364    9771       2279    12.410 %
# 5  112447   66474      33558    29.843 %
# 6  121349  141455      88072    72.577 %
# 7   10096    6276       2962    29.338 %
# 8  124313  135907      90645    72.917 %
# 9   10562    5531       1062    10.055 %

# The network ensemble solves the problem of under-prediction to some extend. 
# Q3 is not too bad. 
# However, the under-prediction of 'G' and 'Hb' is still there. 
# Surprisingly, the buried 'C' is very under-predicted...

# In the 7-state model, the minor classes are about 4%. Network predictions are 
# more or less acceptable. 
# In the 10-state model, the minor classes are about 1~2%. We have very bad 
# network predictions. 
# The HMM filter can solve some problem. But if we have very bad network 
# predictions, the results can not be fixed. 

#######################################
# 7.3 5x2 states network training

# Because the 10-state networks perform badly. We try networks to predict 
# 5-state secondary structures and 2-state residue exposures. 

train_5_state.py
train_2_state.py
# (a1*b1.reshape(2,1)).reshape(1,10)
# Nope. The products are not good enough either. 
# The under-predictions of 'G' and 'Hb' states are still there. 

################################################################################
# 8 YASPIN with 6 states                                                       #
################################################################################

# The 10-states model failed because of the under-predictions of minor states. 
# Try [C E H] x [buried exposed] 6 states. 

# 6 states
#          C  E  H 
# exposed  0  1  2 
# buried   3  4  5

# C 0 661718  0.27088
# E 1 174656  0.07150
# H 2 462873  0.18948
# C 3 341923  0.13997
# E 4 366509  0.15003
# H 5 435196  0.17815

#######################################
# 8.1 network training

train_6_state.py -hu 80 -rs 9

# training error ranges from 1.05940 to 1.24357 in 100 nets 
# the best net is y_6_80_90.net
# the worse one is y_6_80_21.net
test_y6_nets.py 
# y_6_80_90.net    avg_error 1.21127
# ensemble of nets avg_error 1.02938
# ensemble of nets with continuous training targets avg_error 1.03971

#                     single network    network ensemble  continuous targets
#                     avg_out  pos_out  avg_out  pos_out  avg_out pos_out
# Ce 0 661718  0.27088 0.23483 0.44948  0.27871  0.52000  0.26702 0.49252
# Ee 1 174656  0.07150 0.04837 0.21871  0.07573  0.27573  0.07386 0.25695
# He 2 462873  0.18948 0.26609 0.60608  0.19638  0.50854  0.18682 0.47029
# Cb 3 341923  0.13997 0.09763 0.22172  0.12858  0.27679  0.14128 0.28188
# Eb 4 366509  0.15003 0.09800 0.40892  0.14637  0.51291  0.14759 0.50619
# Hb 5 435196  0.17815 0.25508 0.62620  0.17423  0.52481  0.18343 0.52066
res/6_state_nets.png

# with the single network, we have trouble predicting Ee, Ce and Eb classes. 
# with the network ensemble, the false-positives were vastly reduced. 
# However, the Ee and Cb classes are still not optimal. 
# Try train_6_state.py with continous targets
# It produced more or less identical results.
# However, the exposure prediction is better...
# mean squired differences: 0.08862 vs 0.08730
# relative entropy errors:  0.54970 vs 0.54382
train_6_state.py

#######################################
# 8.2 get transition table 

get_6_state_trans_table.py
# the first order tables
# forward transition table
# [0.578967504,0.045589905,0.056470653,0.236287014,0.052066441,0.030618482,],\
# [0.126222160,0.316982537,0.003057262,0.108681372,0.435113044,0.009943625,],\
# [0.077500656,0.000250978,0.491594586,0.048004852,0.002112862,0.380536066,],\
# [0.465856422,0.054522108,0.050676748,0.274507545,0.106307467,0.048129711,],\
# [0.091374957,0.194695969,0.002891795,0.107198818,0.595904023,0.007934437,],\
# [0.043739559,0.000767998,0.415113427,0.032578082,0.005034943,0.502765990,],\

# backward transition table 
# [0.579075021,0.034643395,0.056407173,0.247260689,0.052670994,0.029942728,],\
# [0.166130206,0.316971632,0.000666197,0.105432351,0.408883531,0.001916082,],\
# [0.077599430,0.001152602,0.491576517,0.036954296,0.002290084,0.390427072,],\
# [0.445157025,0.056186849,0.065812502,0.274441631,0.116393401,0.042008591,],\
# [0.090346642,0.207186913,0.002668029,0.097890127,0.595928428,0.005979862,],\
# [0.044732002,0.003986197,0.404555179,0.037313635,0.006680513,0.502732474,],\

test_yaspin_6_esb.py
# Q3 77.71 %
# Q6 61.17 %
# Q2 76.35 % (exposure prediction)

# using the emmission probabilities alone
# Q3 76.47 %
# Q6 58.87 %
# Q2 75.61 %

# not too bad ...

# confidence (marginal likelihood) ranges from 0.2 to 1.00000
# confidence         Q3      Q6      Q2
# 0.00000 - 0.51783  55.89%  33.73%  59.80%
# 0.51783 - 0.60193  70.03%  42.60%  61.70%
# 0.60193 - 0.68635  72.14%  47.44%  65.65%
# 0.68635 - 0.76675  74.54%  52.44%  69.57%
# 0.76675 - 0.83780  76.86%  58.13%  74.29%
# 0.83780 - 0.89489  79.41%  63.24%  77.94%
# 0.89489 - 0.93645  82.45%  69.21%  82.26%
# 0.93645 - 0.96403  85.13%  75.23%  86.81%
# 0.96403 - 0.98211  87.87%  81.02%  90.78%
# 0.98211 - 1.00000  92.73%  88.64%  94.69%

# So far the best is still the old 7-state model
# The 10-state model suffers from the under-prediction of minor classes. 
# The 6-state model losts in secondary structure details. 
# However, to get the exposure predictions and apply higher order Markov model,
# I still want to try the 6-state model. 

# third order transition tables of HMM
third_ord_6_state_trans_table.py

# 216 x 6 transition table 
# 34 of the 216 combinations never observed.
# the largest value is 0.978947, Eb Ee He -> He 

test_yaspin_3order_6state_esb.py

# Q3 77.76 %
# Q6 61.10 %
# Q2 76.31 %
# well, nothing much here either

################################################################################
# summary                                                                      #
################################################################################

#                               Q3         Q2        QX (X can be 6,7 or 10)
# YASPIN 7-state single net     77.78 %              70.96 %
#        7-state net ensemble   78.42 %              71.50 %
#       10-state single net     74.72 %              56.91 %
#       10-state net ensemble   77.94 %              59.79 %
#        6-state net ensemble   77.71 %    76.35 %   61.17 %
#        6-state emmission only 76.47 %    75.61 %   58.87 %
#        6-state 3rd order      77.76 %    76.31 %   61.10 %

# The snow-ball approach of network training failed. 
# Model with more than 7 states/classes failed, because the networks under-
# predict minor classes. 
# Network ensembles work, always better than single networks.
# For the secondary structure prediction, use the old 7-state model. 
# For the residue exposure prediction, use the 6-state first order model. 

################################################################################
# the end                                                                      #
################################################################################
