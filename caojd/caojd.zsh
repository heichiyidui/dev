################################################################################
#                                                                              #
#                              C.A.O.J.D                                       #
#                                                                              #
################################################################################

################################################################################
# Scoring matrices for YASPIN and residue-residue contact prediction           #
################################################################################


################################################################################
# 1. select H class representative domains                                     #
################################################################################
mkdir index

# try get as many alignments as possible for each H class
wc ../cath/bl_out/* | sed 's/..\/cath\/bl_out\///' | grep -v total | \
    awk '{print $1/2 -1 ,$4}' | sort -r -g > t.in

sel_dom_repr.py > index/dom.ls
# 2099 domains
# 4428785 alignments
# 2110 Blast hits per domain.
# The number is smaller than the average number of hits of S35 domains.
# H representatives of the small domain classes tend to have less Blast hits.

# Need to cut alignment number?
# Did some check.
# It seems the class 2 and 4 domains can have deep alignments as well.

################################################################################
# 2. VTML matrix estimation                                                    #
################################################################################

# The method is from T Muller and M Vingron 2000.
# We need the Q (rate matrix) and the Pi matrix (amino acid frequencies)
# Infact, Q = S x Lambda x inv(S), Lambda being the diagonal eigenvalue matrix
# P(t) = S x exp(t x Lambda) x inv(S), for any time t.
# Lambda eigen values are all negative in this case.

# for a alignment with a N matrix, Nij for number of aligned i and j Amino acids
# F for the diagonal Pi for the distribution of amino acids
# argmax SUM_i,j N_i,j log ([ F x exp(tQ) ]_i,j)
# will give us the right t

#######################################
# 2.1 get amino acid frequencies
aa_counter.py
# 296824 residues from 2099 domains

# A 24190  0.08150
# R 15904  0.05358
# N 12453  0.04195
# D 17025  0.05736
# C  3877  0.01306
# Q 11253  0.03791
# E 21524  0.07251
# G 20983  0.07069
# H  6601  0.02224
# I 17336  0.05840
# L 28225  0.09510
# K 17929  0.06040
# M  6204  0.02090
# F 11614  0.03913
# P 13285  0.04476
# S 17158  0.05781
# T 16025  0.05399
# W  3916  0.01319
# Y 10016  0.03374
# V 21306  0.07178

# copy the vtml directory
# Change the Pi file. And it will stay the same this time.

#######################################
# 2.2 alignment distance estimation

# 15 seconds for 3135 alignments.
# about 200 alignment per second

# get_vtml_dis.py

# To run it on nc2, we had to change the script from python3 to python2.7,
# because of the libraries numpy and scipy are not available for python3!

split -l 21 index/dom.ls

for ifile in x?? ; do
    get_vtml_dis.py $ifile > $ifile.out &
done

# should be done in ~ 0.5 hours
cat x??.out > dis.out
# the first round distance estimation:
# min:       0.3118
# mean:      105.9467
# max:       240.8510
# std:       33.8776

################################################################################




################################################################################
# 1. CAO contact clustering

# To reduce the dimension of our problems, we need to further classify the
# residue contacts into groups?


