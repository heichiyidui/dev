################################################################################
#                                                                              #
#                             C.A.T.H                                          #
#                                                                              #
################################################################################

################################################################################
# CATH is a classification of protein domain 3D structures by Orengo CA et al. #
#                                                                              #
# I need CATH for my protein structure prediction projects.                    #
# I need:                                                                      #
#   1. non-redundant domains                                                   #
#   2. residue PSIBlast PSSM                                                   #
#   3. residue secondary structure definitions                                 #
#   4. residue exposure                                                        #
#   5. residue contact definitions                                             #
################################################################################

################################################################################
# 1. Get the CATH files                                                        #
################################################################################

########################################
# 1.1 get the pdb and index files

# the pdb files, 	01-Jul-2013 13:55 	376M
wget http://release.cathdb.info/v4.1.0/CathDomainPdb.S35.tgz
tar xvzf CathDomainPdb.S35.tgz
rm CathDomainPdb.S35.tgz
# 21155 domain files in the dompdb directory

# the index file
mkdir index
wget http://release.cathdb.info/v4.1.0/CathDomainList.S35
mv CathDomainList.S35 index/

awk '{print $1}' index/CathDomainList.S35 | sort  > t1.ls
ls dompdb | sort > t2.ls
diff t1.ls t2.ls
# all domains have corresponding CATH classification

awk '{print $2 "." $3 "." $4 "." $5}' index/CathDomainList.S35| sort | uniq | wc
# 2737 H classes

########################################
# 1.2 remove the transmembrane domains
# The page http://www.cathdb.info/sfam/membrane/ was long gone...
# Searching keyword 'membrane' on CATH gives 348 H families. Way too much.

# Search 'membrane' from the CATH domain description file.
# Searching 'transmembrane' won't pick enough.

wget http://release.cathdb.info/v4.1.0/CathDomainDescriptionFile
grab_description.py | grep membrane | awk '{print $1}' > t1.ls
rm CathDomainDescriptionFile
wc t1.ls
# 308999 domains, 5692 domains with 'membrane' in NAME or SOURCE.
# 316 of them in CATH s35

grab -f t1.ls index/CathDomainList.S35 | awk '{print "rm dompdb/" $1}' > t.out
source t.out
grab -f t1.ls -v index/CathDomainList.S35  > t.out
mv t.out index/CathDomainList.S35
# remove the 316 domains

# 20839 domains left

################################################################################
# 2. QC the CATH files                                                         #
################################################################################

########################################
# 2.1 check residue names are of the standard types
# UNK, PCA, ASX etc residues are to be fixed.

chk_unk.py

# 9 domains
# use psiblast for them.
# 3e2oA02 seq lptXysl -> lptDysl
# 3ze9B00 SEC -> CYS
# 3tguV00 too small, too many UNK, too many breaks, remove
# 3jxvA02 VLKEX -> VLKEG
# 1nthA00 MGVO  -> MGVK
# 2yhxA02 too many UNK, remove
# 1kqfA03 SEC -> CYS
# 2iv2X03 SEC -> CYS
# 2yhxA03 too many UNK, remove


# 4 pdb domain files need to be fixed ...
# using psiblast for it.
# 3e2oA02 seq lptXysl -> lptDysl
# 3jxvA02 seq lkeXegy -> lkeGegy
# 3ed7A00 seq LXXXXPPHG -> -----PPHG
# 2yhxA03 too many UNK, remove

# 20836 domains left

########################################
# 2.2 check atoms with negative occupancy

chk_occu.py

# nothing found this time

########################################
# 2.3 multiple locations of the same atom

chk_multiloc.py > t.in

# 1304 domains with multiple locations of the same atoms

# remove 218 domains with more than 10 atoms with multiloc
# 20618 domains left

chk_multiloc.py > t.in

# 1086 domains need to be fixed now.

# download the original PDB files in case...
mkdir rcsb
awk '{print "wget http://www.rcsb.org/pdb/files/" substr($1,1,4) ".pdb"}' t.in \
    > t.out
sort t.out | uniq > t.sh
# 880 pdb files to be downloaded from RCSB
cd rcsb
source ../t.sh

for ifile in *.pdb; do
    grep "^ATOM " $ifile > t.out; mv t.out $ifile ;
done

cd ..

# for atoms with multiloc atoms, the smallest altCode wins.
mkdir pdb2
rm_multiloc.py
mv pdb2/* dompdb/
rm -r pdb2

########################################
# 2.4 domain size

# check domain sizes
rm t2.out
for ifile in dompdb/* ; do
    echo $ifile >> t2.out ;
    grep " CA " $ifile | wc  >> t2.out ;
done

# then remove 748 domains with > 500 or < 50 residues
# 19870 domains left

########################################
# 2.5 check the atoms are not too close to each other

ls dompdb > t.ls
chk_atom_dis.py t.ls > t.in
# it takes some 80 hours. Use the C++ version

g++ -O4 chk_atom_dis.cpp
a.out t.ls > t.in
# it takes about 0.8 hours.

# 27 files to fix.
# removed 4 domains: 2es7A00 1gjjA00 3zxaC01 2h41A00

########################################
# 2.6 check the alternative residue names with the same chain number id

ls dompdb > t.ls
chk_res_name.py t.ls
# nothing found

########################################
# 2.7 check atoms in the same residue should be close to each other

chk_res_atom_dis.py t.ls

# slow. used cluster or the C++ version
g++ -O4 chk_res_atom_dis.cpp
a.out t.ls > t.out
# Mostly are HH22 atoms popping warning messages.

########################################
# 2.8 check residues (missing CA, N or C atoms)

mkdir pdb2
chkCANC.py t.ls > t.out

# 270 of 19866 domains have missing CA, N or C residues

# remove the two domains with too many breaks
# 1epaA00 160 9 0 9 9
# 1ml9A00 260 10 0 10 10

mv pdb2/* dompdb/
rm -r pdb2

# 19864 domains left

########################################
# 2.9 check for distances between CA and N, CA and C within the same residue

ls dompdb > t.ls
chkResBrk.py t.ls > t.out
# distances should be less than 2.5A and more than 1.0A

awk '{print $2 "_" $3 "_" $4 "_" $5 }' index/CathDomainList.S35 \
    | sort | uniq | wc
# 2515 H classes so far
# 19864 domains

########################################
# 2.10 check for segments and breaks of the mainchain.

#  The distance between the C atom and the N atom of the next residue
#  is mostly (>99%) between 0 and 2.5A.
#  I use 2.5A as the threshold here.

ls dompdb > t.ls
chk_seg_size.py t.ls > t.out
awk '{print $3}' t.out | sort -g | uniq -c

# 10699 domains are single-segment-domains.
# 4140 domains have 1 break.
# 1844 have 2.
# 3181 domains have 3 or more.

# 2.18 segments per domain on average
# mean segment length 70, median 50

# remove domains with segments of less than 20 residues
# 14776 domains left.
# removed 5088 domains.

awk '{print  $2 "_" $3 "_" $4 "_" $5 }'  index/CathDomainList.S35 \
    | sort | uniq -c  | wc
# 2191 H classes

########################################
# 2.11 check the altloc of C, N and CA atoms of the same residue

ls dompdb > t.ls
chkCNCAaltloc.py t.ls > t.out

# Want them to be of the same altloc.
# 25 domains to be fixed.

########################################
# 2.12 mapping residue number

# The residue numbers are messy in PDB files. Different insertion codes,
# chain id changing and alt_loc are troubles to programs like STRIDE and SAP.
# Give them the uniformed numbering without insert or chain id or altloc.
#

ls dompdb > t.ls
res_id_map.py t.ls > index/cath_s35.res_num_map
# save a copy of the residue id mapping

mkdir pdb2

chgResId.py t.ls
rm -r dompdb
mv pdb2 dompdb
# 14776 domains, 2191 H classes left.

# The QC is almost done here.
# But not yet.

################################################################################
# 3. Contact definition                                                        #
################################################################################

#######################################
# 3.1 get residue side-chain and backbone centers.

# For residues with all side-chain heavy atoms, just use the weight center.
# For residues with many missing side-chain heavy atoms, use the projection
# from the C, CA and N atoms.

mkdir ccbc
ls dompdb > t.ls
g++ getCb.cpp -O4
a.out

# 2219103 Cb centers calculated
# 18746 guessed. (less than 1%)

#######################################
# 3.2 get the Delaunay tetrahedralization contact definitions.

# The Delaunay tetrahedralization was performed using TetGen v1.4.
# Use only contacts between side-chain centers, which are not blocked by
# the main-chain to main-chain or side-chain to main-chain contacts.
# Set the distance threshold to 8 A.

mkdir conDef
g++ -c predicates.cxx -O2
g++ -c tetgen.cxx -O2
g++ getTet.cpp predicates.o tetgen.o -O2

ls ccbc > t.ls
a.out

parse_condef.py t.ls > index/cath_s35.condef

rm -r conDef ccbc

#######################################
# 3.3 collect domain sequences
ls dompdb  > t.ls
writeSeq.py > index/cath_s35.seq

#######################################
# 3.4 contact number statistics

contact_number_counts.py

# pretty linear between contact number and domain size
# contact_number = -59.7 + 3.4207 * domain_size
# correlation coefficient is 0.995

# 2237849 residues, 6772808 contacts
# considering one contact involves two residues,
# 6.053 contacts per residue

# contact number per residue:
# max 17, min 0, median 6, mean 6.053
# stddv 2.2753, pretty normal

# residue   avg_number_of_contacts
# 'G'       4.563524
# 'K'       4.908375
# 'E'       5.102626
# 'D'       5.227688
# 'N'       5.282160
# 'R'       5.292399
# 'P'       5.347199
# 'S'       5.377228
# 'Q'       5.385786
# 'T'       5.826501
# 'H'       5.881748
# 'A'       5.889850
# 'C'       6.724250
# 'Y'       6.997959
# 'M'       7.087242
# 'V'       7.107942
# 'W'       7.606361
# 'I'       7.723948
# 'L'       7.740654
# 'F'       7.771821

# It seems more related with hydrophobicity than size.
# T-test and chi-squired test say the different residues have different
# number-of-contacts distributions.

simul_contact_number.py
# We collected the number of contacts each residue could have.
# Now, given a sequence we can randomly draw a set of contact numbers, from the
# distributions we collected. Then a sum contact number can be acquired.
# For the same sequence we run the summing 50 times. With the 50 sum contact
# numbers we can calculate a mean and a standard deviation. The observed domain
# sum contact number gives a Z score, calculated according to this mean and
# standard deviation combination.
#

# The Z score of domains tends to be normally distributed, with a mean of -0.56
# and a standard deviation of 3.04

# It is linearly correlated to domain length.
# Z = -4.8251 + 0.028153 * domain_length
# The correlation coefficient is 0.723.

# After some manual check, any domain with Z score less than -6 will be regarded
# too 'loose' and removed.
# 355 domains to be removed. They are often small domains.

# 14421 domains left.

# update file:

# index/CathDomainList.S35
# index/cath_s35.condef
# index/cath_s35.seq
# index/cath_s35.res_num_map

################################################################################
# 4. DSSP and STRIDE secondary structure definitions                           #
################################################################################

#######################################
# 4.1 run STRIDE and DSSP

# to install STRIDE
# wget http://webclu.bio.wzw.tum.de/stride/stride.tar.gz
# tar xvzf stride.tar.gz
# make
# mv stride ~/bin
#
# and clean the files

# to install DSSP
# wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-i386
# mv dssp-2.0.4-linux-i386 ~/bin/dssp

mkdir stride_ss
awk '{print $1}' index/CathDomainList.S35  > t.ls
awk '{print "stride dompdb/" $1 " > stride_ss/" $1}' t.ls > t.out
source ./t.out

mkdir dssp_ss
awk '{print "dssp dompdb/" $1 " > dssp_ss/" $1}' t.ls > t2.out
source ./t2.out

#######################################
# 4.2 Parse the DSSP and STRIDE definitions

# 89 residues in 87 domains were not assigned any secondary structures by DSSP.
# (double missing assigments in 2 domains.)
# Just treat as their secondary structure are 'X' and their ACC are 'NA'.
# STRIDE has no such problem.

parse_dssp.py > index/cath_s35.dssp
# after change a line, get ACC
parse_dssp.py > index/cath_s35_dssp.acc


parse_stride.py  > index/cath_s35.stride
# or change a line, get ACC
parse_stride.py  > index/cath_s35_stride.acc

rm -r dssp_ss
rm -r stride_ss

#######################################
# to get 3-states definitions:
# H and G can be considered as helix (H), E and B(b) as strand (E)
# and all others as coil (C).

#######################################
# to get 7-states definitions

# Neither STRIDE or DSSP want to assign secondary structures to the begin/end of
# a domain. Most are coil, a few begins (5.9%) are assigned to Helix according
# to STRIDE.

# From 3 states to 7 states:
# C -> C
# H -> H
# E -> E
# CHH -> Hb -> G
# EHH -> Hb -> G
# HHC -> He -> I
# HHE -> He -> I
# HEE -> Eb -> D
# CEE -> Eb -> D
# EEC -> Ee -> F
# EEH -> Ee -> F


################################################################################
# 5. blast the sequences                                                       #
################################################################################

# We need to blast the domain sequences to get PSSMs and multiple sequence
# alignments.

#######################################
# 5.1 download the NCBI nr and the Uniprot uniref 90 databases

nrformat.zsh





#######################################
# 3.1 get the initial sequences
mkdir seq
ls dompdb > t.ls


#######################################
# 3.2 download the nr database

# Downloaded the nr database from ncbi (ftp://ftp.ncbi.nlm.nih.gov/blast/db/).
# Date: 05/10/12
# It was then filtered to remove various non-globular/biased regions using
# the program pfilt (v1.4) by David T. Jones.

cd ../nr
nrformat.zsh
# had some warnings from pfilt
# WARNING - description line truncated - increase BUFLEN!
# so... change BUFLEN from 1048576 (2**20) to 2097152 (2**21) in pfilt.c

cd ../cath

#######################################
# 3.3 get PSSMs

mkdir pssm

# try run blast locally?
ls seq > t.ls
awk '{print "~/bin/psiblast -db ../nr/nr -num_iterations 2 -query seq/" $1\
     " -out_ascii_pssm pssm/" $1}' t.ls > t.out
# 10 mins a job, about 50 or 100 days to finish 14776 jobs on my computer.

init_blast.py t.ls
# 100 sequences done in 8.5 hours. Not too bad. But that's only the first round.
# To get PSSMs, at least two iterations.
# Possible to run on the cluster. But I don't have the cluster access now.

################################################################################
# now we are in trouble...                                                     #
################################################################################

# try remote ncbi blast
mkdir bl_out
ls seq > t.ls
awk '{print "~/bin/blastp -db nr -query seq/" $1 " -out bl_out/" $1 \
    " -max_target_seqs 5000 " \
    " -outfmt \"7 qstart qend qseq sseq \" -remote"}' t.ls > t.sh

# 20 mins per sequence

# Better to try batch searching. 5 sequences in a query file.
# Too many sequences leads to SIGXCPU (24) error. NCBI wants jobs to be finished
# in less than one hour (?) combined CPU time.

ls bl_in > t.ls
awk '{print "~/bin/blastp -db nr -query bl_in/" $1 " -out bl_out/" $1 \
    " -max_target_seqs 5000 " \
    " -outfmt \"7 qstart qend qseq sseq \" -remote"}' t.ls > t.sh

# submit 5 jobs simultaneously.
# 1.1 mins per sequences. About 1300 sequences per day.
# Need 12 Days to finish 14776 sequences.

# Run psiblast and get PSSM from NCBI is problematic.
# The " -num_iterations " option is not compatible with "-remote "
# The PSSM from NCBI is different and got '-I -I ... -I' lines.
# It's likely that the PSSM direct from NCBI is using the query sequence only.

################################################################################

# All sequences have hits, often around 5000.
# We asked for '-max_target_seqs 5000', but sometime we got more than 5000 hits.
# In these cases, the hits are often very duplicated.
# All aligned segments have capital letters. In addition to 'B', 'Z' and 'X',
# some are 'J', 'U' and 'O'.

# De-duplication. Stack the pairwise alignments together.
# Replaced leading and ending '-'s with '.'s.
# Top 5000 hits if more are found.

mkdir bl_out2

parse_blast_out.py

# 14776 domains blasted. 2237849 residues.
# The average residue reading depth is 2311.
# 365 residues reading depths are zero.











init_blast.py t.ls
# sometime got

split -l 100 t.ls
# put the commandline into t.sh
# init_blast.py $1

# t.sh:
#-----------------------------
# #!/bin/sh
# #$-S /bin/sh
# #$ -cwd
# #$ -q bignode.q,long.q,short.q
# #################################################
#
# id_ls=$1
# ./init_blast.py $id_ls
#-----------------------------

for lsfile in x?? ; do
    qsub -l h_vmem=10G t.sh $lsfile
done

# It takes 8.5 hours to blast 100 sequences

mkdir cseq
calign.py t.ls
# was using kalign
# got error messages in kalign such as
# 25504 Bus error  ~/bin/kalign c_align_in/1avvA00 > c_align/1avvA00
# Failed about 100 times out of 16,556 domains.
# Used MUSCLE instead.
# MUSCLE is faster, but the output order is not stable.

# anyway, got the c-seq

# for comparisons, get the domain sequences with segments seperated with 'x'
mkdir bseq
ls dompdb > t.ls

writeBseq.py t.ls

# 16556 domains, 1-26 breaks per domain
# 9357 have no break, 3310 have 1, 1471 have 2, etc.
# domain length 14-1146, average 148, median 130

# 34624 segments, length 1-779, average 71, median 50
# 18068 breaks

######
# from cseq,  1-25 breaks per domain, average 0.94
# 9690 no break, 3344 1, 1502 2 etc.
# 32144 segments, length 1-759, average 76, median 57
# Of the 16556 domain sequences, 15832 were modified with added residues.
# That's 95.6%.

compBreaks.py t.ls
# with default 0.5 cuts
# (The query sequence has a gap at the point, and the most occuring amino acid
# is more than half in the column of the alignment)

# Of the 18068 real domain breaks,
# 6365 have not been assigned (35%)
# 11703 were assgined via calign.py
# calign.py found 15588 gaps to be filled.
# 11703 of them are corresponding to real breaks in domains (75%).
# a lot of them are 'm' (5088) and 'x' (2360)
rm -r bseq

#######################################
# 3.4 get PSSM with the second round of Blast

# use mostly default parameters.
# 3 iterations

ls cseq > t.ls
awk '{print "~/bin/psiblast -db ../nr/nr -num_iterations 3 -query cseq/" $1\
     " -out_ascii_pssm pssm/" $1}' t.ls > t.out
split -l 100 t.out

# header ...
head -n 4 t.sh > t.header
for ifile in x?? ; do
    cat t.header $ifile > t.out;
    mv t.out $ifile
done

for psifile in x?? ; do
    qsub -l h_vmem=10G $psifile
done

getPssm.py t.ls
# or
split -l 100 t.ls
for lsfile in x?? ; do
    qsub -l h_vmem=6G t.sh $lsfile
done

# 100 blast jobs took some 18-22 hours

# two short cseq (1g3iW02 and 1kyiS02) got not blast hits
# Using ncbi blast web server
# It automatically adjust parameters for short input sequences to get some hits.
# So, add options of our psiblast searches for the two
# -evalue 200000 \
# -word_size 2 \
# -matrix PAM30 \
# -gapopen 9 -gapextend 1 \
# -comp_based_stats 0 \
# -inclusion_ethresh 0.005 \

# 2696068 PSSM for 2696068 residues in cseq
ls pssm > t.ls
parse_pssm.py > index/cath_s35.pssm
rm -r pssm

# 53921360 values
# min -16, max 13, mean -1.86, sdv 2.98, median -2
# sounds like -2 +- 3


################################################################################
# Done                                                                         #
################################################################################
