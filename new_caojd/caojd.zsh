################################################################################
#                              C.A.O.J.D                                       #
################################################################################

################################################################################
#                                                                              #
# CAOJD is a method to predict protein residue contacts from mulitple sequence #
# alignments.                                                                  #
#                                                                              #
# We need stricter domain QC.                                                  #
# We need the CAO and YASPIN results.                                          #
# I decided to remove the gap-filling procedure.                               #
#                                                                              #
# What information to use from a multiple sequence alignment?                  #
# 1. PSSM                                                                      #
# 2. predicted secondary structures and contact numbers                        #
# 3. point residue conservation                                                #
# 4. adjusted correlated covariation between residue pairs                     #
# 5. contact patterns between secondary structure segments                     #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
# 1. Extra QC of CATH domains                                                  #
# 2. YASPIN preparation                                                        #
# 3. CAO and VTML estimation                                                   #
# 3. YASPIN network training                                                   #
# 4. CAOJD                                                                     #
#                                                                              #
################################################################################

################################################################################
# 1. Extra domain QC                                                           #
################################################################################

################################################################################
# 1.1 domain breaks 
cp ../cath/index/CathDomainList.S35 index/
awk '{print $1}' index/CathDomainList.S35 > t.ls
chk_chain_brk.py t.ls > t.in 

# we want domains with 4 or less segments, and 60 < domain length < 500
awk '{if ($2<500 && $2>60 && $3 < 5) print $1}' t.in > t.ls 

grab -f t.ls  index/CathDomainList.S35 > t.out
mv t.out index/CathDomainList.S35 
# 2870 domain removed 

grab -f t.ls t.in > t.out 
mv t.out t.in 

# remove domains with short segments with less than 16 residues
awk '{for (i=4;i<=NF;i++) if ($i<16) print $1}' t.in | sort | uniq > t.ls
# another 1924 domains to remove 
grab -f t.ls -v t.in > t.out
mv t.out t.in 

grab -f t.ls -v index/CathDomainList.S35 > t.out
mv t.out index/CathDomainList.S35 

# 11762 domains left

################################################################################
# 1.2 DSSP vs STRIDE

awk '{print $1}' index/CathDomainList.S35 > t.ls

mkdir stride_ss
awk '{print "~/bin/stride ../cath/dompdb/" $1 " > stride_ss/" $1}' t.ls > t.out

t.head 
#-------------------------------------- 
#!/bin/sh
#$-S /bin/sh
#$ -N stride
#$ -cwd
#$ -q long.q,bignode.q,short.q
#--------------------------------------
# end of t.head 

split -l 100 t.out 
for ifile in x?? ; do 
    cat t.head $ifile > t.out 
    mv t.out $ifile 
    qsub $ifile 
done 
# done in seconds 

mkdir dssp_ss
awk '{print "~/bin/dssp -i ../cath/dompdb/" $1 " > dssp_ss/" $1}' t.ls > t.out
split -l 100 t.out 
for ifile in x?? ; do 
    cat t.head $ifile > t.out 
    mv t.out $ifile 
    qsub $ifile 
done
# done in seconds

# put domain sequences into index/cath_s35.seq in FASTA two lines format
# 1803014 AA left

awk '{print $1}' index/CathDomainList.S35 > t.ls
parse_stride.py t.ls 
# index/cath_s35.stride index/cath_s35.stride.acc

parse_dssp.py t.ls 
# index/cath_s35.dssp   index/cath_s35.dssp.acc

# 76 domains have residues missing assignments in DSSP definitions. 
# These residues are only the C terminal residue with missing "O" atoms. 
# Use STRIDE definitions in such cases. 

#   STRIDE            DSSP
# B  17362        B  18452
# b    659        b      0
# C 326127        C 371318
# E 407501        E 391938
# G  67445        G  65936
# H 626855        H 601618
# I    275        I    321
# S      0        S 157009
# T 363775        T 203407

# The difference between STRIDE and DSSP secondary structure definitions 
# can be quite big, even when we treat all DSSP 'S' as 'T'

# remove the 150 domains with difference larger than 30% 
# 11612 domains left in index/CathDomainList.S35

aac_coef.py 

# Remove 1351 domains with less than 0.99 correlation coefficient between 
# DSSP acc and STRIDE acc. 
# 10261 domains left. 

# linear regression of sum_acc vs domain_length shows
# sum_acc ~= 2749.6 + 37.155 * domain_length 

# remove 52 domains with residue larger than 3879.5
# now residue  -4 < Z < 4
# 10209 domains left.
# sum_acc ~= 2786.6 + 36.758 * domain_length

################################################################################
# 1.3 contact numbers

awk '{print $1}' index/CathDomainList.S35 > t.ls
cont_count.py > t.in 

# use ../cath/index/cath_s35.condef 
# the number of local contacts can be very linear vs the number of residues
# beta residues should have 2 each
# alpha residues should have 4 (or 5?) each

res/loc_cnt_vs_resi_number.agr

# each residue should have about 1.6 x 2 or more global contacts
# smaller domains might have less

res/glo_cnt_vs_resi_number.agr

# now remove the (many) domains with less that 1.6 global contacts per residue
# remove 3972 domains

# 6237 of 16556 domains are left. 62%, 10319 domains are gone!
# 992 of 2679 H families are left.
# That's a lot removed!

################################################################################
# QC ends here                                                                 #
################################################################################

################################################################################
# 2. YASPIN preparation                                                        #
################################################################################

################################################################################
# 2.1 secondary structure and exposure states transition tables

# 6237, 1187589 residues

# STRIDE definitions:
# B   13255  0.01116
# E  285404  0.24032
# G   47000  0.03958
# H  352654  0.29695
# C  245445  0.20668
# I     216  0.00018
# T  243615  0.20513

# B, E -> E
# G, H -> H
# C, I, T -> C 

# C  489276  0.41199
# E  298659  0.25148 
# H  399654  0.33653

# All domains start and end with 'C'.

# 'CHH','EHH' -> 'G'
# 'HHC','HHE' -> 'I'
# 'HEE','CEE' -> 'D'
# 'EEH','EEC' -> 'F'

# C  489276  0.41199
# D   53061  0.04468
# E  192537  0.16212
# F   53061  0.04468
# G   44744  0.03768
# H  310166  0.26117
# I   44744  0.03768

# saved into index/cath_s35.stride.ss7

# SS7_to_i={'C':0,'D':1,'E':2,'F':3,'G':4,'H':5,'I':6}

# 1181352 transitions between states

#######################################
# 2.1.1 secondary structure transitions 
# transition table:
#  C        D        E        F        G        H        I
# [0.784431,0.106882,0.022797,0.000000,0.085883,0.000006,0.000000,],
# [0.000000,0.000000,0.871751,0.128249,0.000000,0.000000,0.000000,],
# [0.055885,0.000000,0.699642,0.240245,0.004228,0.000000,0.000000,],
# [0.953921,0.000000,0.000000,0.000000,0.046079,0.000000,0.000000,],
# [0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,0.000000,],
# [0.000010,0.000000,0.000000,0.000000,0.000000,0.855732,0.144258,],
# [0.955413,0.032027,0.012560,0.000000,0.000000,0.000000,0.000000,],

# reverse transition table:
#  C        D        E        F        G        H        I
# [0.784431,0.000000,0.022276,0.104787,0.000000,0.000006,0.088500,],
# [0.972993,0.000000,0.000000,0.000000,0.000000,0.000000,0.027007,],
# [0.057194,0.240245,0.699642,0.000000,0.000000,0.000000,0.002919,],
# [0.000000,0.128249,0.871751,0.000000,0.000000,0.000000,0.000000,],
# [0.927163,0.000000,0.018192,0.054644,0.000000,0.000000,0.000000,],
# [0.000010,0.000000,0.000000,0.000000,0.144258,0.855732,0.000000,],
# [0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,0.000000,],

#######################################
# 2.1.2 residue exposure transitions

# max residue acc Values from Sander & Rost, (1994), Proteins, 20:216-226
# MAX_ACC['A']=106.0; MAX_ACC['R']=248.0; MAX_ACC['N']=157.0; MAX_ACC['D']=163.0
# MAX_ACC['C']=135.0; MAX_ACC['Q']=198.0; MAX_ACC['E']=194.0; MAX_ACC['G']= 84.0
# MAX_ACC['H']=184.0; MAX_ACC['I']=169.0; MAX_ACC['L']=164.0; MAX_ACC['K']=205.0
# MAX_ACC['M']=188.0; MAX_ACC['F']=197.0; MAX_ACC['P']=136.0; MAX_ACC['S']=130.0
# MAX_ACC['T']=142.0; MAX_ACC['W']=227.0; MAX_ACC['Y']=222.0; MAX_ACC['V']=142.0

# RSA (relative solvent accessible area) ranges from 0 to 1.8
# median 0.251
# mean 0.305 std 0.277
# 11.2 % are 0

index/cath_s35.stride.rac

# use 0.25 as the threshold for RSA
# rsa > 0.25 is exposed, 50.1 %
# rsa < 0.25 is buried,  49.9 %
# This time, after much stricter QC on domains, RSA is smaller. 
# The threshold 0.25 is generally accepted. 

# together with 3-states secondary structure 
# C buried  0 179246  0.15093      
# E buried  1 207654  0.17485      
# H buried  2 206464  0.17385      
# C exposed 3 310030  0.26106
# E exposed 4  91005  0.07663
# H exposed 5 193190  0.16267

index/cath_s35.stride.ss6

# forward transition table
# [0.289975, 0.111109, 0.048042, 0.450053, 0.055359, 0.045462],
# [0.108339, 0.613232, 0.007219, 0.086004, 0.182438, 0.002769],
# [0.035512, 0.005943, 0.523777, 0.044269, 0.000901, 0.389598],
# [0.252407, 0.059374, 0.030777, 0.558178, 0.048637, 0.050627],
# [0.110016, 0.446250, 0.009901, 0.120960, 0.309752, 0.003121],
# [0.054521, 0.002759, 0.403566, 0.081510, 0.000254, 0.457389],

# backwrad transition table
# [0.290019, 0.125665, 0.040956, 0.428599, 0.055926, 0.058836],
# [0.095804, 0.613232, 0.005909, 0.086919, 0.195571, 0.002567],
# [0.041663, 0.007260, 0.523777, 0.045315, 0.004364, 0.377620],
# [0.265058, 0.058744, 0.030064, 0.558129, 0.036209, 0.051797],
# [0.108917, 0.416285, 0.002044, 0.162464, 0.309752, 0.000538],
# [0.042135, 0.002976, 0.416367, 0.079663, 0.001470, 0.457389],

################################################################################
# 2.2 PSSM from PSIBlast search

# using 10/Jan/2015 NCBI NR database

mkdir pssm 

~/bin/psiblast -db ../nr/nr -num_iterations 3 -query seq/in.seq \
-out_ascii_pssm pssm/out.pssm 
# 100 sequences took some 50 hours on the cluster.

parse_pssm.py

index/cath_s35.pssm

################################################################################
# 
# Try redo Jens' structural alphabet. 

# get the triplets CA distances
# see chk_chain_brk.py

# mean CA-CA distance 3.80356 std 0.05291
# why the heck we have 2.5 and 4.5?

# The 2.5 distances are bad. They are just there in the pdb domain files. 
# Someone didn't fix their PDB files properly...

# removed the distances beyond [3.6,4.0]
# 10915 removed out of 1159437 triplets

# mean 3.8061, std 0.0251
# normalize it. 

# The VBEM package from M Beal is not properly working on this distances. 
# It finishes clustering the 3D dots with 8 or 7 components, with the largest 
# component occupies some 70%. That's way too big. I'm expecting at least one 
# component for a secondary structure group. 

# The 3D plots of the distance triplets show nothing... 

# library(Bmix)
# library(mvtnorm)
# y <- as.matrix(read.table('xaa',header=FALSE))
# alpha <- 4; rho <- 0.8
# params <- c(c(0,0,0), #gamma
#               0.1, #kappa
#               3, #nu
#               3, #gam0
#               diag(.5,3) #psi0
#               )
# res <- mix(y, alpha=alpha, g0params=params, times=NULL, rho=rho, cat=0,
#            state=NULL, read=FALSE, print=FALSE, N=11486, niter=0)
# The clustering using Bmix is another disappointment.
# The package is too slow. 
# So far, the MCMC procedure is suggesting 20 clusters. 

# library(bgmm)
# y <- as.matrix(read.table('xaa',header=FALSE));
# res=unsupervised(X=y,k=13);
# write.table(res$mu,file='clu_13.out',\
#     append=TRUE,row.names=FALSE,col.names=FALSE)
# 
# Well, the three distances are correlated. 
# Cut the input into 100 pieces, the 13 x 100 centers are grouped around 5 or 
# more clusters on the diagnal. 
# Will not see that if use 20 centers in the algorithm.

# Use the python scikit learning tools instead. 
# With the full data, VBGMM and DPGMM are both giving non-sense result. 
# Use the simple GMM with diagnal (or full?) covariances instead. 
# BIC suggests some 11 (?) clusters. 
gmm.py 

# got 16 clusters with diagnal covariances
gmm_res 

# or 10 clusters with full covariances

# Overall, trying to rebuild the structural alphabet failed.
# The three distances between Ca atoms alone will not be sufficient.  

################################################################################
# 3. CAO and VTML estimation                                                   #
################################################################################


################################################################################
# 4 CAO estimation                                                             #
################################################################################

# After QC, we have 2024 H families left.
# Most have only a few members, 30 H families have 50 or more members. 
# For each of the 2024 domains, in case the number of pairwise alignments is 
# very big, we randomly select about 1000 of the alignments to avoid bias. 

################################################################################
# 4.1 select domains

sort -g -k 2 -r index/cath_s35.aln_depth | awk '{print $1}' > t.ls
# want domains with deeper alignments be there first
sort_table -f t.ls index/CathDomainList.S35 | \
    awk '{print $1, $2 "." $3 "." $4 "." $5}'  > t.in

# in python 
#--------------------------------------
#!/usr/bin/env python3 

ifile=open('t.in')
h_classes=[]
for line in ifile:
    cols=line.split()
    if h_classes.count(cols[1]) > 1:
        continue
    h_classes.append(cols[1])
    print(cols[0])
ifile.close()

#--------------------------------------

index/cath_dom_rep.ls
# End up with 3123 domains in 2024 H families. Each family has 1 or 2 domains. 

################################################################################
# 4.2 VTML estimation

# method from T Muller and M Vingron 2000
# We need the Q (rate matrix) and the Pi matrix (amino acid frequencies)
# Infact, Q = S x Lambda x inv(S), Lambda being the diagonal eigenvalue matrix 
# P(t) = S x exp(t x Lambda) x inv(S), for any time t.
# Lambda eigen values are all negative in this case. 

# for a alignment with a N matrix, Nij for number of aligned i and j Amino acids
# F for the diagonal Pi for the distribution of amino acids
# argmax SUM_i,j N_i,j log ([ F x exp(tQ) ]_i,j)
# will give us the right t

#######################################
# 4.2.1 VTML alignments transfer to matrices
# sum alignments into number matrices

mkdir vtml_aln 
get_vtlm_aln.py index/cath_dom_rep.ls
# 3657320  alignments
# average 1171 alignments per domain 

#######################################
# 4.2.2 VTML distance estimation 

# copy the old CATH VTML matrices into vtml/

mkdir vtml_dis 
/share/apps/R_3.0.2/bin/Rscript vtml_dis.R index/cath_dom_rep.ls
# 100 domains each, 30 jobs 
# about 3 hours 

# 3657320 distances
# min 1.0  max 225
# mean 102

# 4156 alignments with distance "1.00000002058295". We need to ignore them.
# Ignore any alignments with distances less than 1.5.

#######################################
# 4.2.3 stack alignment matrices together

sum_vtml_aln.py > vtml_aln.in
rm -r vtml_aln 
rm -r vtml_dis 

#######################################
# 4.2.4 estimate R and Q and Pi 

# get the distances of the sum matrices
#/share/apps/R_3.0.2/bin/Rscript on the cluster
Rscript vtml_sum_dis.R 

# estimate Q and Pi
Rscript get_vtml_pi_q.R

cp Q2 vtml/Q 
cp Pi2 vtml/Pi 

# Repeat many iterations till the Pi values converge.

# Pi changed, especially on tryptophan.
# The estimated frequency of tryptophan rised from 1.25 to 4.95%. 
# Dawn J. Brooks et al. 2002 esitmated amino acid frequencies in the Last 
# Universal Ancestor (LUA) according to PAM. The result is very similar to what
# we got here. 

#######################################
# 4.2.5 write new VTML matrices
Rscript write_vtml.R 

mv vtml50 vtml75 vtml100 vtml150 vtml/

# This differencies between the new and old VTML 100 are in the range of
# -0.40 ~ 0.40 , mainly on residue D and C.

################################################################################
# 4.3 CAO estimation

#######################################
# 4.3.1 CAO alignments transfer to matrices

mkdir cao_aln 
get_cao_aln.py t.ls 

# The jobs run reasonably fast. But they created 1.4T of files!

#######################################
# 4.3.2 CAO distance estimation 

/share/apps/R_3.0.2/bin/Rscript cao_dis.R index/cath_dom_rep.ls
# This is going to take about 3 mins per alignment 
# 21 years. 

# Need to stack the CAO alignments together according to VTML distances. 
# It is going to be less accurate. But I don't have 21 years for this. 

# To stack multiple alignments:
# update 't_dis_split.ls' according to new 'vtml_sum.dis'

# to estimate vtml distance of an alignment matrix xaa.vmat 
/share/apps/R_3.0.2/bin/Rscript vmat_dis.R xaa.vmat

# stack a small collection of domain alignments
stack_cao_aln.py xaa

random_shuffle_lines.py index/cath_dom_rep.ls > t.ls
split -l 60 t.ls
for ifile in x?? ; do qsub t.sh $ifile ; done 

# stacking them use about 8 hours with 53 jobs on the cluster 
# one list failed. Didn't bother about the 60 domains

# sum them together to 100 x 160000 matrix in cao_aln.in 

# 1527993410 / 2 aligned contact pairs observed. 
# 4775 counts per element in a 160000 matrix

# two 0 entries found in the sum of the 100 matrices 
# CW <-> HP never observed. 

# the largest entry is 84210 with 14081732 counts
# LL <-> LL is the most frequently observed entry. 
# then LV <-> LV
# then VL <-> VL
# then IL <-> IL 

#######################################
# 4.3.3 CAO estimation

# get the distances of the sum matrices
#/share/apps/R_3.0.2/bin/Rscript on the cluster
Rscript cao_sum_dis.R 

# estimate Q and Pi
Rscript get_cao_pi_q.R

cp Q2 cao/Q 
cp Pi2 cao/Pi 
# repeat... several iterations till converge

#######################################
# 4.3.4 write new CAO matrices
Rscript write_cao.R 
mv cao100 cao150 cao200 cao250 cao/

# the difference between VTML + VTML and CAO is
# C-C is way more conservitive than what we expected from VTML.
# W-K (or any W contacts) are less conservitive than what we expected from VTML.

# 168 is the CAO distance of the 50th alignment matrice. 
# The correlation coefficient is 0.84 
# between the observed cao168 pmat and the estimated cao168 pmat

################################################################################
# 5. contact pattern definition                                                #
################################################################################

#######################################
# 5.1 plot the contact maps

plot_cont_mat.py 
