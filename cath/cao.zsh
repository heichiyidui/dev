################################################################################
#                                                                              #
#                               C.A.O                                          #
#                                                                              #
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
# The difference is in CAO we are looking at residue side-chain contacts.      #
#                                                                              #
################################################################################

################################################################################
# 1. multiple sequence alignment of BLAST hits                                 #
################################################################################

#######################################
# 1.1 blast consensus sequences 

ls cseq > t.ls 
blast_cseq.py t.ls 
# subject range should be at least half of query length

# use the cluster for it. 
vi t.sh 
#--------------------------
# #!/bin/sh
# #$-S /bin/sh
# #$ -cwd
# #$ -q bignode.q,long.q,short.q
# #######################################
# 
# id_ls=$1
# ./blast_cseq.py $id_ls
#--------------------------

split -l 100 t.ls 
for lsfile in x?? ; do
    qsub -l h_vmem=10G t.sh $lsfile
done

# around 6 hours to blast 100 sequences

ls cseq > t.ls 
awk '{print "~/bin/muscle -quiet -in c_align_in/" $1 " > c_align/" $1}' t.ls \
 > t.out 
mkdir c_align

split -l 100 t.out 
head -n 4 t.sh > t.header
for ifile in x??; do cat t.header $ifile > t.out; mv t.out $ifile; done 
for ifile in x??; do qsub -l h_vmem=10G $ifile; done 

#######
# aligning 2-9047 sequences
# min 2, max 9047, mean 366, sdv 317, median 384
# MUSCLE uses around 10 hours to align 9047 sequences of avg length 102.

ls c_align > t.ls 
mkdir seq_aln 

parse_maln.py

################################################################################
# 2. SAP alignments                                                            #
################################################################################
# randomly select domain pairs from the same H families 
# of the 2679 H classes, 1252 have 1 domain, 460 have 2, 242 3 etc.
# 

select_sap_pairs.py > t.sh 
# from the 1427 H classes, selected 17076 pairs of domains 
# Not too bad, the 3 most popular of the 11785 domains were selected 10 times. 
# 3467 domains selected once
# 2397 domains selected twice
# 1875 trice
# 1624 4
# 1288 5
#  728 6
#  310 7
#   75 8
#   18 9
#    3 10
# not too biased towards certain H classes. 

mkdir sap_aln 
source t.sh 

# or to use the cluster 
vi t.header
#!/bin/sh
#$-S /bin/sh
#$ -cwd
#$ -q bignode.q,long.q,short.q

split -l 1000 t.sh 
for ifile in x?? ; do  cat t.header $ifile > t.out; mv t.out $ifile; done 
for ifile in x??; do qsub $ifile ; done 

################################################################################
# 3. parse SAP alignments                                                      #
################################################################################

ls sap_aln > t.ls 

parse_sap_aln.py 
# 17076 SAP pairwise alignments
# 506 with alignment scores less than 200 and will be ignored.

# 2123693 aligned residue pairs, in the 16570 alignments.
# 421855 with pairwise score less that 1.0  (20%)
# 501437 with pairwise score less that 2.0  (24%)
# use the 1622256 residue pairs with >= 2.0 scores

# 6099170 alignments
# 16570 * 2 are SAP alignments and 6066030 are sequence alignments
# 0.54% is SAP, that's too few? 

################################################################################
# don't bother re-estimate VTML and CAO                                        #
# The old estimation with CATH 3.5 should be good enough for quite a while.    #
################################################################################

################################################################################
# 4. VTML estimation                                                           #
################################################################################

# 4.1 estimate alignment distances 
# 4.2 calculate R, Q and Pi 
# 4.3 go back to 4.1 

#######################################
# 4.1 alignment distance estimation 
# method from T Muller and M Vingron 2000
# we need the Q (rate matrix) 
# infact, Q = S x Lambda x inv(S), Lambda being the diagonal eigenvalue matrix 
# P(t) = S x exp(t x Lambda) x inv(S), for any time t.
# Lambda eigen values are all negative in this case. 

# for a alignment with a N matrix, Nij for number of aligned i and j Amino acids
# F for the diagonal Pi for the distribution of amino acids
# argmax SUM_i,j N_i,j log ([ F x exp(tQ) ]_i,j)
# will give us t

ls seq_aln > t.ls
mkdir n_mat

getAaNmat.py

mkdir cao_dis
#R CMD BATCH cao_dis.R
R  --vanilla < cao_dis.R t.ls  
# around 10 hours to estimate distances 
# split into 25 jobs, each deal with around 100 domains
# much faster on the cluster 

#######################################
# 4.2 calculate R and Q and Pi 

# 4.2.1 sum n_mats into 100 bins according to cao distances 
ls cao_dis > t.ls
sum_n_mat.py > t.mat 

# 4.2.2 get dis of the 100 sum n_mats 
R CMD BATCH cao_s_dis.R

# 4.2.3 get R, Q and Pi, with alpha set to 0.02 
R CMD BATCH getR.R
mv Q2 Q
mv Pi2 Pi 

#######################################
# 4.3 save the Q, Pi and VTML50 VTML75 and VTML100 
mkdir vtml 

R CMD BATCH  writeVTML.R  
cp Q Pi vtml/
cp vtml100 vtml75 vtml50 vtml/

# the order of the amino acids is 
#'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'

################################################################################
# 5. CAO estimation                                                            #
################################################################################

#######################################
# 5.1 get CAO pi matrix 

ls seq_aln > t.ls 
getCaoPi.py > Pi 

# 505397 global contacts in 2429 domains, 320225 residues
# average 1.58 contact per residue, a residue is in 3.16 contacts 
# global means sequence distance larger than 4 
# say, 1-6, 1-7 etc 

plotCaoPi.R 
# it turns out C-C con_tacts are way more than what we expect 
# larger and hydrophobic amino acids are more in contacts 
# for charged amino acids, basic avoid basic, acid avoid acid 
# there are more basic-acid contacts. 
# say for D, there are more D-R, D-H and D-K, less D-E and D-D 

#######################################
# 5.2 get the 400 x 400 alignment N_mat matrices 

sum_cao_n_mat.py  > t.mat 
# about 4 hours to sum them up 
 
# some types of contact mutations are never observed ,
# W-C, W-W, W-M, C-C, W-H 

#######################################
# 5.3 get CAO R, Q and Pi 
# for the moment, use the old VTML distances

# add Laplacian prior
R CMD BATCH  getCaoR.R 
cp Q2 Q
cp Pi2 Pi 

#######################################
# 5.4 update distances 

# add Laplacian prior
R  --vanilla < sum_cao_dis.R

#######################################
# 5.5 re-estimate Q and Pi 

R CMD BATCH getCaoR.R 

# save the results
mkdir  cao 
mv Q2  cao/Q 
mv Pi2 cao/Pi 

################################################################################
# 6. Clustering of contact patterns                                            #
################################################################################

################################################################################
# It seems the clustering of binary contact patterns is necessary. I assume    #
# predicting a contact between two residues is possible. However, predicting   # 
# contacts between two groups of adjacent residues could to be more            #
# practical. And the predictions should be more informative.                   #
################################################################################

########################################
# 6.1 preparing the sub-matrices to be clustered. 

# add secondary structures into it. 
# 3 types of secondary structures, C, H and E
# 9 combinations. 
# 10 million ExE ? 

# overall, we have around 263,840,050 sub_matrices to cluster
# most of which are all 0 

# use H class representatives, remove all zero matrices
# use DSSP secondary structures
# 5x5 matrices: fxf
# 4,778,559 none-zero matrices

get_sub_matrices.py > fxf.out
shuffle_sub_mat.py fxf.out > fxf.in ; rm fxf.out 

# 889561 C C 
# 807795 H H
# 655896 C H
# 600141 H C
# 480553 C E
# 458289 E C
# 375673 E E
# 265518 E H
# 245133 H E

# of the 4,778,559 5x5 non-zero sub-matrices
# number of contacts ranges between 1 and 13 (observed once)
# there are 136880 types (way less than 2**25 - 1)

# 27,184,803 zero sub-matrices
# 5469007 C C 
# 4527696 H H 
# 4125849 C H 
# 4060742 H C 
# 2467295 C E 
# 2476780 E C 
# 1591014 E E 
# 1271164 E H 
# 1195256 H E 

# got all 25  patterns with 1 contact,  occurances 81K ~ 109K
# got all 300 patterns with 2 contacts, occurances 400 ~ 16K
# got 2295  of  2300 patterns with 3 contacts. (occurances 1 ~ 9K)
# got 10791 of 12650 patterns with 4 contacts. (occurances 1 ~ 6K)

# 164 patterns with contact number 11, mostly once, 3 patterns observed twice
# max contact number is 13 (observed once), 12 (18 patterns, once each)

sed 's/[CHE]  [CHE]  //' fxf.in > t.in; mv t.in fxf.in 

head -n 25000 fxf.in > t.in 

#######################################
# try the Bayesian Hierarchical Clustering methrod 
# the bclust package
install.packages('bclust')

# on the cluster
/share/apps/R_3.0.2/bin/R

install.packages('bclust',lib='~/R') 
library('bclust',lib.loc='~/R') # to load in on the cluster

init_cluster.R
# the results are rather disappointing 
# 1. it's very slow
# 2. it often cluster mostly unrelated patterns together
# 3. it produced way too many clusters. 

#######################################
# try the Bayesian k-means 

wget https://sites.google.com/site/kenichikurihara/academic-software/\
bayesian-k-means/bkm_src.tar.gz
tar xvzf bkm_src.tar.gz  
cd bkm_src 

DATA=load('t.in');
DATA=DATA';

[z,k] = bkm_bu(DATA) 
# still producing way too many clusters. 
# still too slow.
# shall I put some noise into the data? 

#######################################
# variational Bayesian mixture factor analysis 

wget http://www.gatsby.ucl.ac.uk/~beal/software/vbmfa/vbmfa.tar.gz
tar xvzf vbmfa.tar.gz 
cd vbmfa 

matlab 
data=load('t.in')';
[data ppp] = preprocess(data);
net=vbmfa(data,1,0,0,10); 
# 27 components given 500  3 x 3 sub-matrices
# 78 components given 2000 3 x 3 sub-matrices
# way too many 

#######################################
# k-means

split -l 50000 fxf.in 
Rscript km.R xaa 
# should put the centers of 50-means clustering into xaa.out 

# for ifile in x??; do  Rscript km.R $ifile; done 

t.sh 
#--------------------------
#!/bin/sh
#$-S /bin/sh
#$ -cwd
#$ -q bignode.q,long.q,short.q
#################################################
Rscript km.R $1
#---------------------------

for ifile in x??; do qsub t.sh $ifile; done 
cat x??.out > fxf.cents

# supper fast

#######################################
# try clues to find the number of cluster

install.packages('clues',lib='~/R') # to install on the cluster

library(clues,lib.loc='~/R') # to load on the cluster
library(clues)

cx <- as.matrix(read.table('fxf.cents'))
b <- clues(cx,strengthMethod='CH')
write.table(file='fxf_clues.clu',b$mem,row.names=FALSE,col.names=FALSE)

# 12 clusters on 3 x 3 sub-matrices according to clues CH index
# 24 clusters on 5 x 5 sub-matrices, 
# 49 clusters on 7 x 7 , 
# not right... 

#######################################
# try BIC 
library(mclust)
cx <- as.matrix(read.table('fxf.cents'))
fit2 <- Mclust(cx,G=1:100)
summary(fit2)

# within cluster variance is too small
# way too many clusters 

########################################
# 6.2 hierarchical clustering 

head -n 25000 fxf.in > t.in
# in R 
cx <- as.matrix(read.table('t.in'))
d <- dist(cx, method = "manhattan") # distance matrix
fit <- hclust(d, method="ward")

plot(fit, labels=FALSE) # display dendogram
rect.hclust(fit, border="red", k = 66 ) 
groups <- cutree(fit, k= 66 )

write.table(file='t.clu',groups,row.names=FALSE,col.names=FALSE)

# very fast
# to cluster 25,000 matrices, use about 7G  memory
# to cluster 50,000 matrices, use about 30G
# to cluster 60,000 matrices, use about 40G

# it seems hierarchical clustering is the only reasonable one. 
# I choose the number of clusters to be 66.

# on the cluster
# use the new R
/share/apps/R_3.0.2/bin/
# It is way faster than using R_2.14.1

# tried iterative clustering of all non-zero sub-matrices with HC
clu_def.py > t2.pat ; grep -v dis t2.pat > t.pat 
clu_def.py > t3.pat ; grep -v dis t3.pat > t.pat 
# ...
# the clustering will not converge 
# maybe because of using the Manhattan distance
# It will jump from A to B, then back to A. 

# apply k-means again with the HC patterns, 
# on the whole set of 4.8 million sub-matrices

x   <- as.matrix(read.table('fxf.in'))
pat <- as.matrix(read.table('t.pat'))

kc <- kmeans(x,clu,iter.max = 100)
write.table(file='pat_66.txt',kc$centers,row.names=FALSE,col.names=FALSE)

# 14k to 153k sub-matrices per patterns
# average about 100k

# 66 patterns, plus an all-zero one. 
# we end up with 67 patterns.
index/pat_67.txt

########################################
# 6.3 plot the 5x5 patterns

# put 5x5 0 into t00.dat

get_mat_patterns.py 

# in matlab 

clims = [ 0 1 ];
colormap(1-gray);
subplot(7,10,1 ); t=load('t00.dat') ; imagesc(t,clims)
subplot(7,10,2 ); t=load('t01.dat') ; imagesc(t,clims)
subplot(7,10,3 ); t=load('t02.dat') ; imagesc(t,clims)
# ...
subplot(7,10,67); t=load('t66.dat') ; imagesc(t,clims)

mkdir cont_patterns
# save the plot into cont_patterns/pat_67.png

########################################
# show some complete contact maps:

>1t6oA00
XXXXXXXXCCCHHHHHHHHHHSCSCHHHHHHHHHHHHHCCSHHHHHHHHHHHHHHHCXXX
XXXXXXXXCCCHHHHHHHHHHCCCCHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHCXXX
# up-down-up helix , parrell and anti parellize

>1s3rA02
XXXXXXXXXCEECCCCCCCEEEEEEEETTEEEEEEEEEEEEEEECCCEEEEEEEEEEEEEEC
XXXXXXXXXCEECCCCCCCEEEEEEEETTEEEEEEEEEEEEEEETTTEEEEEEEEEEEEEEC
# 3 strands up-down-up with parellal and anti-parella

>2j7jA01
XXXXXXXXCEECCSTTCCCEESSHHHHHHHHHHHHCXXXXXXXXX
XXXXXXXXCEEECTTTTCEEECCHHHHHHHHHHHHCXXXXXXXXX
# 2 strands and 1 helix

# The DSSP and STRID definitions are almost identical.

# plot the contact maps to cont_patterns/

################################################################################
# The End                                                                      #
################################################################################
