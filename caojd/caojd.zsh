################################################################################
#                                                                              #
#                              C.A.O.J.D                                       #
#                                                                              #
################################################################################

################################################################################
# Scoring matrices for YASPIN and residue-residue contact prediction           #
################################################################################

################################################################################
#                                                                              #
# CAO (contact accepted mutation).                                             #
# Kuang Lin, Jens Kleinjung, William R Taylor, Jaap Heringa 2003               #
# Computational Biology and Chemistry  2003, vol 27 issue 2 Pages 93-102       #
#                                                                              #
# It's after the PAM (Point Accepted Mutation), a Markov model of amino acid   #
# replacements in proteins.                                                    #
#                                                                              #
# The VTML method is from T Muller and M Vingron 2000.                         #
# The difference is in CAO we are looking at residue side-chain contacts.      #
#                                                                              #
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
# because the numpy and scipy libraries were not available for python3!

random_shuffle_lines.py  index/dom.ls > t.ls
split -l 21 t.ls

for ifile in x?? ; do
    get_vtml_dis.py $ifile > $ifile.out &
done

# should be done in ~ 0.5 hours
cat x??.out > dis.out
rm x??
rm x??.out

# to check the distribution of distances:
grep -v ">" dis.out > t.ls
listdis t.ls > t.dat

# the first round distance estimation:
# min:       0.3118
# max:       240.8510
# mean:      105.9467
# std:       33.8776

# ...

# the eighth round
# min:       0.3091
# max:       252.5416
# mean:      108.2719
# std:       35.2503

# the ninth round
# min:       0.3091
# max:       252.5146
# mean:      108.2730
# std:       35.2509

#######################################
# 2.3 sum up the alignment matrices, get the new Q

mkdir sum_mat
get_vtml_Q.py
# Don't forget to change the MAX_DIS constant!
# It takes about a hour.

# and then go back to distance estimation.
# C-E and C-K might have negative values in the Q matrix when alpha is small.

# Stop here after the ninth round.

#######################################
# 2.4 cluserting of amino acids:
# using -log(pt120), in R

# library('fastcluster')
#
# pt120=as.matrix(read.table('sum_mat/050',header=F))
# dis=as.dist(-log(pt120),diag=F,upper=F)
#
# Labels=c('A','R','N','D','C','Q','E','G','H','I',
#          'L','K','M','F','P','S','T','W','Y','V')
#
# tree = hclust(dis, method='ward.D2' )
# plot(tree,labels=Labels)
#
# clu = rect.hclust(tree,k=5)

#######################################
# 2.5 get VTML scoring matrices

get_vtml.py 100
# 100 here is the PAM distance
cp t.out vtml/vtml100

# To plot it in R:
#
# a=as.matrix(read.table('t.out'))
# Labels=c('A','R','N','D','C','Q','E','G','H','I',
#          'L','K','M','F','P','S','T','W','Y','V')
# rownames(a)=Labels
# colnames(a)=Labels
# heatmap(a)

################################################################################
# 3. CAO matrix estimation                                                     #
################################################################################

#######################################
# 3.1 CAO contact frequency

# 886435 contacts
# According to the amino acid pair frequencies, the sequence distance clearly
# matters.

# Sequence distance groups: 1, 2, 3, 4, 5 and more
# or 1 and 2, 3 and 4, 5 and more ?

mkdir cao
get_cao_pi.py

# Save the Pi matrices into cao/Pi_1, cao/Pi_2, cao/Pi_3, cao/Pi_4 and cao/Pi_5.

# to plot the log-odds of frequency matrices in R

# Labels=c('A','R','N','D','C','Q','E','G','H','I',
#          'L','K','M','F','P','S','T','W','Y','V')
#
# freq_aa = as.matrix(read.table('vtml/Pi'),header=F)
# base_aa_aa = freq_aa %*% t(freq_aa)
# rownames(base_aa_aa)=Labels
# colnames(base_aa_aa)=Labels
# heatmap(base_aa_aa,symm=T)
#
# freq_aa_aa = as.matrix(read.table('cao/Pi_5'),header=F)
# rownames(freq_aa_aa)=Labels
# colnames(freq_aa_aa)=Labels
#
# freq_aa_aa = log(freq_aa_aa) - log(base_aa_aa)
# heatmap(freq_aa_aa,Colv = "Rowv")
#
# # write.table(freq_aa_aa,file='t.out', col.names=F,row.names=F)

# remember the frequency matrices are not symmetric.
# For Pi_5, there are basically two groups: big and hydrophobic, and others.
# For Pi_4, G and P (maybe C) are special.
# For Pi_3, C is special.
# For Pi_2, unclear, rather 3 groups.
# For Pi_1, two groups, grouping is similar to Pi_5. The frequencies are very
# different.

#######################################
# 3.2 sum up cao alignment matrices

# using vtml distances to sum CAO alignments into 100 matrices
get_cao_sum_mat.py
# Note the matrices 94, 96, 97 and 98 have no alignment data.
# we used the outer product of pi instead.

#######################################
# 3.3 CAO Q esitmation

get_cao_Q.py
# a few minutes

get_cao_dis.py > cao.dis
# 30 to 50 minutes

# the first round cao distances:
# quartiles: 1.7617 93.1874 202.8065 318.2812 1897.4218
# the second round
# quartiles: 1.7677 95.6950 213.6681 343.6129 2043.7113
# the third round
# quartiles: 1.7699 95.9406 214.8316 346.4199 2089.7011

# after the first round, calculate the Q and the distances using only the first
# 85 sum alignment matrices. The last summary matrices are of too few samples.

# The problem is still under-sampling.
# Too many bins to fill.

#######################################
# 3.4 CAO estimation

get_cao.py 200

# now t.out is the CAO200 matrix
# min:       -4.5226
# max:       6.6471
# mean:      -1.0756
# std:       1.3066
# quartiles: -4.5226 -1.9866 -1.2318 -0.3605 6.6471

# better use VTML100 and CAO200 for alignments

# From 3669406 pair-wise alignment of PAM distance 50-150,
# we obtained the frequencies of aligned contact-pairs.
# From CAO model, we estimated the frequencies of aligned contact-pairs at
# CAO distance 200.

# The log of two sets of frequencies were compared in
pic/cao_vs_obv.pdf
# It is quite linear. The correlation coefficient is 0.96703.

################################################################################
# 4. classify contacts                                                         #
################################################################################

# sum sum_mat/040 ~ sum_mat/060
# now all bins are covered.
# log of the sum_mat is pretty normal

cont_cluster.R

# group contact types into 20 clusters.
# we used hierarchical clustering with manual selected number of cluster.

cont_cluster.py > t.in

# the second part of
cont_cluster.R

# the plots
pic/cont_groups_1.pdf
pic/cont_groups_2.pdf
# show the clustering of contacts.

################################################################################
# The End                                                                      #
################################################################################


