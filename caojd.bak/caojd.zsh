################################################################################
#                                                                              #
#                              C.A.O.J.D                                       #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
# CAOJD is a method to predict protein residue contacts from mulitple sequence #
# alignments.                                                                  #
#                                                                              #
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
# 1. QC of domains                                                             #
################################################################################

#######################################
# 1.1 We use the domain list from YASPIN. 
# It went through some QC of transmemberane removal and secondary structure 
# consistency check. Both useful here. 

~/bin/random_shuffle_lines.py ../yaspin/index/s.ls > index/s.ls

#######################################
# 1.2 get the new contact definition file, keep the global contacts only
get_global_conts.py > index/cath_s35.condef

#######################################
# 1.3 count residue global contact numbers 
count_res_cont_num.py

#######################################
# 1.4 simulate the number of contacts according to the domain length and 
# residue contact number distribution 
simu_dom_cont_num.py

# 518 domains with Z scores less than -7 
# These domains are often long singleton helix. 
# remove these 518 domains, left with 15934 domains in the list 

# removed 44 domains in the second round after updating the frequencies 
# left with 15890 domains

# now 14.5% of residue have no global contacts

#######################################
# 1.5 blastp again 

# We want deeper multiple alignments 
# installed blast-2.2.29+ for faster batch searching 

# use -max_target_seqs 5000 to get at most 5000 gi hits

blastp -db ../nr/nr -outfmt "7 sseqid sstart send slen pident" \
    -max_target_seqs 5000 -query t1.seq > t1.out

# use 15G memory 

# found blastp is over-reporting sequences? 
# over-lapping segments or even the identical segments are reported? 

# lots (8918) "Selenocysteine (U) at position XXX" warning messages.
# If there's 'U' in the sequences, change them to 'C'

# It took some 2~3 hours to finish searching 100 sequences. 
# 2 domains have 0 hits. Manually remove them to avoid troubles in parsing.
# 1472 have less than 200 hits. They are to be removed.
# 6123 have 5000 or more hits. That's about the half. 
# 52527664 segments his altogether. 
# The sequence identity ranges from 15 to 100 %. The mean is 40%. 

# Parse the blast hits files. Use only one segment (the one with the highest 
# sequence identity) for each gi.

parse_bhits.py xaa.out > xaa.in 

mkdir bseqs

# to extract GI 4139528 4-330
blastdbcmd  -db ../nr/nr  -entry 4139528 -range 4-330 -outfmt  \
    "%s" >> bseqs/12asA00

awk '{print "blastdbcmd -db ../nr/nr -entry " $2 "\
 -range " $3 "-" $4 " -outfmt \042%f\042 >> bseqs/" $1}' xaa.in > xaa.bdb 

for ifile in x?? ; do 
    awk '{print "blastdbcmd -db ../nr/nr -entry " $2 "\
 -range " $3 "-" $4 " -outfmt \042%f\042 >> bseqs/" $1}' ${ifile}.in \
        > ${ifile}.bdb 
done 

t.head
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N bdb  
#$ -cwd
#$ -q long.q,bignode.q,short.q
#--------------------------------------

for ifile in x?? ; do
    cat t.head ${ifile}.bdb > t.out ; mv t.out ${ifile}.bdb 
done 

rm bseqs/*
for ifile in x??.bdb ; do qsub $ifile ; done 

# slow 
# try -entry_batch option 
mkdir batch_files 

#######################################
# 1.7 get the bseqs to be aligned
mkdir temp
ls bseqs > t.ls

add_dom_seq.py t.ls 
# cut the headers of fasta hits
# add the domain sequence to the head 

rm -r bseqs
mv temp bseqs
mkdir temp

#######################################
# 1.8 kalign or muscle?

# kalign is about 5 times faster.
# muscle shuffles the sequences around. However, with the domain id there, it
# should be easy to fix it back. 
# muscle is too slow. Go kalign 

ls bseqs | awk '{print "~/bin/kalign bseqs/" $1 " > temp/" $1}' > t.ls
split -l 100 t.ls 
for ifile in x?? ; 
    do cat t.head $ifile > t3.out ; 
    mv t3.out $ifile ; 
done 

for ifile in x?? ; do qsub $ifile ; done 

# 669 sequences failed to align. 
# Might be kalign found invalid letter 'J' in the sequences. 
# 'J' ('XLE') is Leucine or Isoleucine.

# replaced those 'J' with 'L' and retry kalign. Done.

#######################################
# 1.9 remove duplications in kaln

# 'B' will be replaced with 'N'
# 'Z' will be replaced with 'Q'
# 'X' will be replaced with '-'
# 'O' will be replaced with 'K'
# 'U' will be replaced with 'C'
# then aligned sequences duplications will be removed. 

# then termini GAPs will be replaced with 'u'
# then any columns with GAPs in the domains will be removed. 
# remove lines with 90% or more GAPs

mv temp kaln 

mkdir temp 
ls kaln > t.ls

parse_kaln.py t.ls 

# reduced the size from 47587108 to 7265408 
# 85% removed...

wc temp/* | sed 's/temp\///' | grep -v total | sort -g > t.ls 

awk '{if ($1<200) print $0}' t.ls | wc
# 2001 domains with less than 200 sequences in the alignments will be removed...
awk '{if ($1<200) print "rm temp/"$0}' t.ls > t.out 
source t.out 
mv temp seq_aln 
# 13887 domains left

# - 288522179  / 5765504497 = 0.050043 
# u 1071535712 / 5765504497 = 0.185853 
# A 410491974  / 5765504497 = 0.071198 
# C 59427702   / 5765504497 = 0.010307 
# D 264204139  / 5765504497 = 0.045825 
# E 285740919  / 5765504497 = 0.049560 
# G 351459715  / 5765504497 = 0.060959 
# F 169354798  / 5765504497 = 0.029374 
# I 272317888  / 5765504497 = 0.047232 
# H 106869120  / 5765504497 = 0.018536 
# K 219003901  / 5765504497 = 0.037985 
# M 95297152   / 5765504497 = 0.016529 
# L 419656346  / 5765504497 = 0.072787 
# N 165781828  / 5765504497 = 0.028754 
# Q 154119449  / 5765504497 = 0.026731 
# P 189427610  / 5765504497 = 0.032855 
# S 244846982  / 5765504497 = 0.042468 
# R 234508379  / 5765504497 = 0.040674 
# T 236225015  / 5765504497 = 0.040972 
# W 53245938   / 5765504497 = 0.009235 
# V 335983729  / 5765504497 = 0.058275 
# Y 137484022  / 5765504497 = 0.023846 

# 32843423 ncbi nr sequences been aligned to them 
# 2365 sequences per domain

# domain wise mean coverage is 77.7%
# mean sequence identity in covered region is 39.2%

# residue wise some got 0 coverage. average is 75.9%. 

#######################################
# 1.10 transpose sequence alignments 
mkdir temp

transpose_aln.py

rm -r seq_aln 
mv temp seq_aln 

# 13887 domains 
# length ranges from 25 to 1146, mean 156.8 
# align depth ranges from 200 to 4957, mean 2366

#######################################
# 1.11 get PSSM 

# put domain sequences into seq/ 

mkdir pssm 

awk '{print "~/bin/psiblast -db ../nr/nr -num_iterations 3 -query seq/" $1\
     " -out_ascii_pssm pssm/" $1}' index/s.ls > t.out 

split -l 100 t.out 
for ifile in x?? ; do cat t.head $ifile > t.out ; mv t.out $ifile ; done 

for ifile in x?? ; do qsub $ifile ; done
# Checked the 13887 PSSM files. They are fine. 

################################################################################
# 2 check multiple sequence alignments to find signals for contact prediction  #
################################################################################

#######################################
# 2.1 residue conservation 

res_cons_prod.cpp 

# this is the fastest to get of all features. 
# try different gap penalty combinations to seperate contacts vs non-contacts

# The conservation of residues can be measured with the vtml100 matrix 
# all-against-all score. 
# use residue conservation score products

# After normal survival function normalization, the scores should be normal. 
# Remember the non-contacts is the overwhelming majority (~98%).
# The only matter is how far the mean of contacts scores from zero. 

g++ res_cons_prod.cpp -O4 -o res_cons_prod
# -O4 here matters
# needs around 1.6G memory

for gap_vtml in  -1.7 -1.6 -1.5 -1.4 -1.3 ;  do 
    for end_vtml in -0.5 -0.4 -0.3 -0.2 -0.1 -0.05 -0.02 -0.01 ; do 
        res_cons_prod $gap_vtml $end_vtml > r_${gap_vtml}_${end_vtml}.out 
    done 
done 

# 0.418909 in r_-1.3_-0.1.out is the largest one 
res_cons_prod -1.3 -0.1 > r_-1.3_-0.1.out

# Tast time the optimal combination wass -1.5 , -0.3
# This time we got better alignments.

res/res_cons_rank_norm.png
# looks good. However, 
# remember the contacts are relatively few ... only 1.792% of all positions 
res/real_res_cons_rank_norm.png

#######################################
# 2.2 test basic cao score 

test_basic_cao.cpp
# g++ -O4 make it way faster. 

split -l 45 index/h.ls 
# 28 lists of 45 domains 

for gap_cao in -2.5 -2.0 -1.5 -1.0 -0.5 ; do 
    for end_cao in  -0.20 -0.10 -0.05 -0.02 -0.01 ; do 
        for ifile in x?? ; do 
./test_basic_cao $ifile $gap_cao $end_cao  > ${ifile}_${gap_cao}_${end_cao}.out
        done 
    done 
done 

for gap_cao in -5.5 -5.0 -4.5 -4.0 -3.5 -3.0 -2.5 -2.0 ; do 
    for end_cao in  -1.5 -1.3 -1.1 -0.9 -0.7 -0.5 -0.3 ; do 
        for ifile in x?? ; do 
             qsub t.sh $ifile $gap_cao $end_cao  
        done 
    done 
done 

#         -5.5    -5.0    -4.5    -4.0    -3.5    -3.0    -2.5    -2.0
# -1.5    0.29174 0.29144 0.29109 0.29062 0.29030 0.29003 0.29025 0.29044
# -1.3    0.28925 0.28895 0.28853 0.28798 0.28751 0.28716 0.28715 0.28693
# -1.1    0.28675 0.28641 0.28595 0.28538 0.28478 0.28426 0.28406 0.28394
# -0.9    0.28401 0.28370 0.28327 0.28259 0.28212 0.28150 0.28110 0.28078
# -0.7    0.27986 0.27947 0.27915 0.27848 0.27791 0.27704 0.27648 0.27599
# -0.5    0.27509 0.27446 0.27383 0.27332 0.27257 0.27170 0.27061 0.26981
# -0.3    0.26850 0.26776 0.26705 0.26616 0.26520 0.26393 0.26252 0.26111

for gap_cao in -5.5 -5.0 -4.5 -4.0 -3.5 -3.0 -2.5 -2.0 ; do 
    for end_cao in -3.5 -3.3 -3.1 -2.9 -2.7 -2.5 -2.3 -2.1 -1.9 -1.7 -1.5; do 
        for ifile in x?? ; do 
             qsub t.sh $ifile $gap_cao $end_cao  
        done 
    done 
done 

#       -5.5   -5.0    -4.5    -4.0    -3.5    -3.0    -2.5    -2.0
# -3.5  0.30364 0.30407 0.30508 0.30555 0.29040 0.30491 0.30368 0.30228
# -3.3  0.30247 0.30273 0.30293 0.30299 0.30396 0.30436 0.30350 0.30217
# -3.1  0.30172 0.30203 0.30207 0.30203 0.30235 0.30315 0.30287 0.30156
# -2.9  0.30044 0.30058 0.30049 0.30031 0.30046 0.30224 0.30252 0.30112
# -2.7  0.29959 0.29948 0.29944 0.29899 0.29902 0.30012 0.30135 0.30042
# -2.5  0.29825 0.29807 0.29794 0.29783 0.29753 0.29800 0.28502 0.29926
# -2.3  0.29741 0.29734 0.29714 0.29699 0.29676 0.29688 0.29943 0.29873
# -2.1  0.29614 0.29586 0.29577 0.29541 0.29524 0.29508 0.29655 0.29723
# -1.9  0.29490 0.29468 0.29464 0.29427 0.29407 0.29391 0.29483 0.29683
# -1.7  0.29327 0.29294 0.29267 0.29237 0.29204 0.29193 0.29214 0.29322
# -1.5  0.29175 0.29144 0.29110 0.29062 0.29031 0.29003 0.29025 0.29045

# failed to find the best seperation
# use gap_cao -2.5 and end_cao -2.3 for the moment

# the score is not that stable. 
# go back with domain QC. 

# tell the script to print domain id and mean contact score 
split -l 50 index/s.ls 

for gap_cao in  -2.0 ; do 
    for end_cao in  -1.5; do 
        for ifile in x?? ; do 
             qsub t.sh $ifile $gap_cao $end_cao  
        done 
    done 
done 

# most mean CAO scores are between -0.5 and +1.0
# When the score is higher than +1.0, quite often we have an extended domain 
# with only a few global contacts. 
# When the score is lower than -0.5, we can often find large chunks of gaps in 
# the alignments. If not so, the domains are often extended as well. 

# lesson: 
# 1. we need to remove domains with few global contacts. 
# 2. we need to remove alignment with many gaps. 
# 3. we need to remove over-conserved alignments
# 4. we need to remove short domains

#######################################
# 2.3 additional domain QC 

# 13887 domains in index/s.ls 

###################
# 2.3.1 remove domains with less than 60 residues or more than 500 residues
# 601 domains removed. 13286 domains left. 

###################
# 2.3.2 rerun simu_dom_cont_num.py 
count_res_cont_num.py 
simu_dom_cont_num.py 

# remove domains with scores larger than +5.0 or less than -5.0
# 858 domains removed. 12428 domains left in the first round 
# 44 domains removed, 12384 domains left in the second round 
# still 13 domains with scores larger than 5.0 in the third round. 
# Didn't bother

# more or less 4 global contacts per residue 

# nothing special about sequence alignment column
# sequence identity and gap percentage distributions

# alignment
# sequence identity should be 15~70%
# gap percentage should be < 50%
# when more than 40% columns have gap >50%, the domain should be removed. 
# when more than 30% columns have gap+identical AA > 90%, remove it. 

# awk '{if ($3>0.15 && $3<0.7 && $4<0.5 && $5<0.4 && $6<0.3) print $0}' t.in \
# > t2.in

# 10229 domains left, 2155 domains removed. 
# that's 17.4% of alignments gone

# 10299 domains in index/s.ls


#######################################
# 2.4 re-run 2.1 and 2.2 
split -l 45 index/h.ls 
# 33 lists

for gap_cao in -5.5 -5.0 -4.5 -4.0 -3.5 -3.0 -2.5 -2.0 ; do 
    for end_cao in -3.5 -3.3 -3.1 -2.9 -2.7 -2.5 -2.3 -2.1 -1.9 -1.7 -1.5; do 
        for ifile in x?? ; do 
             qsub t.sh $ifile $gap_cao $end_cao  
        done 
    done 
done 
# 2904 jobs
# each about 1 hour or less 


for gap_cao in -5.5 -5.0 -4.5 -4.0 -3.5 -3.0 -2.5 -2.0 ; do 
    for end_cao in -5.5 -5.0 -4.5 -4.0 ; do 
        for ifile in x?? ; do 
             qsub t.sh $ifile $gap_cao $end_cao  
        done 
    done 
done 

# do have small difference between the contacts and non-contacts
# on all-againist-all CAO scores
res/cao_rank_norm.png
#      contact non-contacts
# mean 0.2701  -0.0104 
# std  1.0940   1.0064

################################################################################
###### FUCK! I got the mat_len wrong!                                     ######
################################################################################




#######################################
# 3.3 test permutation cao score

permu_cao_rank.cpp

# randomly permute alignment 100 times, re-calculate CAO scores
# for non-contacts, 51.5 % orginal CAO score is better than all permute scores.
# for contacts, 58.9% 
# for other rankings, non-contacts is always higher. 

# calculate contact z scores according to the permutation. 
permu_cao_z.cpp
# the program takes forever ... 

res/permu_cao_z.png
#      contact  non-contacts
# mean 0.252822 -0.0114706
# std  1.09268   0.995415

#######################################
# 2.4 write domain aln feature file:
dom_aln_feature.cpp 
normal_dom_feature.cpp 

# 62610355 in the test set 
# 1136896 contacts  1.815%
# 129229099 in the train set 
# 2295050 contacts  1.776%

grep -v ">" test.in | \
    awk '{if ($5==1 || ($5==0 && rand()<0.02)) print $0}' > t.out 
mv t.out test.in 

grep -v ">" train.in | \
    awk '{if ($5==1 || ($5==0 && rand()<0.02)) print $0}' > t.out 
mv t.out train.in 

# 1136896/2365961 in test 
# 2295050/4832573 in train 

train_cao_dom_fea.py 

################################################################################
# 3 try plot the first neural network prediction results.                      #
################################################################################

#######################################
# 3.1 plot the contact maps

# realized that we need the domain length in the contact definition matps

res/g_conts_vs_dom_length.png
# on average, 4 global contacts per residue 

# got the out-liners in the domain contact numbers 
simu_dom_cont_num.py






################################################################################
# TODO: residue contact number prediction                                      #
################################################################################

#######################################
# 2.1 fit the residue global contact numbers to two distributions
count_res_cont_num.py | head -n 15 > t.dat 

R
library(mixdist)
res_cont_num_data=as.matrix(read.table('t.dat',header=FALSE))

fit <- mix(res_cont_num_data, mixparam(c(0.5,3,6),.333), 'binom',
    mixconstr( consigma = "BINOM", size = c(14, 14, 14)))

# got 6 warning messages
fit
png('res/res_cont_num_dist.png')
plot(fit)
dev.off()

fit$parameters
#          pi        mu     sigma
# 1 0.1726840 0.3596348 0.5919429
# 2 0.4420848 2.7280904 1.4820543
# 3 0.3852312 5.6165786 1.8339291

python3 
import scipy, scipy.stats
x = scipy.linspace(0,14,15)
pmf_1 = scipy.stats.binom.pmf(x,14,0.0256882)
pmf_2 = scipy.stats.binom.pmf(x,14,0.1948636)
pmf_3 = scipy.stats.binom.pmf(x,14,0.4011842)

pi_1=0.1726840
pi_2=0.4420848
pi_3=0.3852312 

pmf=pmf_1*pi_1 + pmf_2*pi_2 + pmf_3*pi_3

# pretty good fit 

#######################################
# 2.2 residue global contact number prediction
train_res_cont_num_net.py 

# 100 networks 
# training error ranges from 0.85917 to 0.93684
# 15 to 93 iterations

