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
#   1. non-redundant, QCed domains                                             #
#   2. residue PSIBlast PSSM and multiple alignments of them                   #
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

# The expected sum contact number is linearly correlated to the domain length.
# C = -4.8251 + 0.028153 * domain_length
# The correlation coefficient is 0.723.

# The Z scores of real sum contact numbers tend to be normally distributed,
#  with a mean of -0.56 and a standard deviation of 3.04 .

# After some manual check, we decided that domains with Z scores less than -6
# are to be regarded as too 'loose' and removed.
# 355 domains to be removed. They are often small domains.

# 14421 domains left.
# 2204188 residues
# Median domain length 133, mean domain length 153.
# 6695153 contacts, 6.07 contacts per residue

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
# STRIDE has no such problem.

parse_dssp.py > index/cath_s35.dssp
# after change a line, get ACC
parse_dssp.py > index/cath_s35_dssp.acc


parse_stride.py  > index/cath_s35.stride
# or change a line, get ACC
parse_stride.py  > index/cath_s35_stride.acc

# Eventually used STRIDE assignments for those 'X's and 'NA's in DSSP.

rm -r dssp_ss
rm -r stride_ss

# DSSP       STRIDE        DSSP_vs_STRIDE
# H 735171   H 766522      HH 720327
# G 79946    G 81348       GG 66905
# E 476838   E 495802      EE 471225
# B 22740    B 22321       BB 17564
# C 450516   C 395994      CC 288656
# T 248147   T 441881      TT 173017
# S 190393
# I 437      I 320         II 192

# DSSP_ACC                STRIDE_ACC
# min:       0.0000       min:       0.0000
# max:       379.0000     max:       404.1000
# mean:      55.0408      mean:      55.1434
# std:       50.5024      std:       49.9600

#######################################
# to get 3-states definitions:
# H and G can be considered as helix (H), E and B as strand (E)
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

mkdir seq
# Then write domain sequence files there.

mkdir pssm

# to get one PSSM

~/bin/ncbi-blast-2.3.0+/bin/psiblast \
    -db ../nr/uniref90 \
    -num_iterations 2 \
    -save_pssm_after_last_round \
    -query seq/1ul4A01 \
    -out_ascii_pssm pssm/1ul4A01

# do it for 14421 sequences

# 50 jobs, each has 289 sequences. All are finished in 4 days, on NC2.

#######################################
# 5.3 blast for multiple sequence alignments

mkdir bl_out

~/bin/ncbi-blast-2.3.0+/bin/blastp \
    -db ../nr/uniref90 \
    -max_target_seqs 5000 \
    -outfmt "7 qstart qend qseq sseq sseqid"\
    -query seq/1ul4A01 \
    -out bl_out/1ul4A01

# do it for 14421 sequences

# 10 jobs, each has 1443 sequences. Most are finished in 4.5 days, on NC2.


#######################################
# 5.4 via NCBI remote blast

# It can also be done with the NCBI blast server
# ls seq > t.ls
# awk '{print "~/bin/blastp -db nr -query seq/" $1 " -out bl_out/" $1 \
#     " -max_target_seqs 5000 " \
#     " -outfmt \"7 qstart qend qseq sseq sseqid \" -remote"}' t.ls > t.sh

# 20 mins per sequence

# Better to try batch searching. 5 sequences in a query file.
# Too many sequences leads to SIGXCPU (24) error. NCBI wants jobs to be finished
# in less than one hour (?) combined CPU time.

# Can submit 5 jobs simultaneously.
# It takes 1.1 mins per sequences, about 1300 sequences per day.
# We need 12 Days to finish 14776 sequences.

# To run psiblast and get PSSM from NCBI is problematic.
# The " -num_iterations " option is not compatible with "-remote "
# The PSSM from NCBI is different and got '-I -I ... -I' lines.
# It's likely that the PSSM from NCBI is built using the query sequence only.

# With the version 2.3.0+, PSIBlast can now store the PSSM after the first
# iteration. The PSSM might be more useful now.

#######################################
# 5.5 parsing PSSM files

parse_pssm.py > index/cath_s35.pssm
# gzip it to save some space

# PSSM values:
# min:       -16.0000
# max:        13.0000
# mean:       -2.1046
# std:         3.2161
# median and mode:     -2

# The histogram is roughly normal, with a bit positive skewness.
# Can be normalized like (x + 2) / 4

#######################################
# 5.6 parsing blast output

# All sequences have hits, often around 5000.
# We asked for '-max_target_seqs 5000', but sometime we got more than 5000 hits.
# In these cases, there are multiple hits from the same sequences.

# All aligned segments have capital letters. In addition to 'B', 'Z' and 'X',
# some are 'J', 'U' and 'O'.

# De-duplication (use the first hit of the same sequence only).
# Stack the pairwise alignments together.
# Replace non-standard residues with '-'.
# Only include subject sequences which are less than 100% identical.
#
# The median of the number of hits is 2960.

mkdir bl_out2

parse_blast_out.py

# 14421 domains, 2204188 residues
# 38248057 subject sequences have been aligned to the domains.
# After parsing, the median per domain is 2544.
# The mean is 2652.
# reading depth of residue:
# min 0, max 5000, mean 2311, median 1931, std 1816

# alignment coverage: (the portion of domain been aligned to)
# min 0.038, max 1    , mean 0.7924, median 0.845
# alignment percent identity: (in the aligned part)
# min 0.111, max 0.998, mean 0.3936, median 0.3697, std 0.108

rm -r bl_out
mv bl_out2 bl_out

################################################################################
# 6. SAP structure alignments                                                  #
################################################################################

# We need some domain structure alignments for the clustering of residue-residue
# contacts or local structure.

select_sap_pairs.py

# 2099 different CATH H classes
# 1221 H classes with two or more domains
# selected 2957 pairs of domains to be aligned.

mkdir sap_aln
# ~/bin/sap dompdb/3lbeB00 dompdb/3r87A00 > sap_aln/3lbeB00_3r87A00
# ...

parse_sap_aln.py > index/cath_s35.sap_aln
# 2923 pairs past
# on average, there are 141 residues in the smaller domains,
# and 100 aligned residue pairs in each alignment

# the mean sequence identity is 20%
# the medean is 17.7%.

# 22 alignments have 100% identity. But that shouldn't matter.

################################################################################
# Done                                                                         #
################################################################################
