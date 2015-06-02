################################################################################
# QC of the batch 3 expression data                                            #
################################################################################

################################################################################
# 1. prepare input files                                                       #
################################################################################

cd /home/kuang/dev/batch_3_expression

#######################################
# 1.1 the sample table

# expression sample ids from Genome\ studio/BATCH3\ 91014/sampletable.txt 
# 564 ids in the expression data

# phenotype data sample ids from 
# Subject_Demographics_April2015_FINAL.ods

# In the pheno table ids, somehow '0's were removed. 
# Like '9370786090_A' was changed to '9377869_A'. 
# Fortunately, it is a one-to-one mapping. Fixed the pheno table ids. 
# One id '94244436_H' was not found in the expression sample table. 
# I guessed the real id is '9402444036_H'.

# 531 rows in the pheno table

# All ids in the expression and pheno table match now.

# 48 chip ids, 46 have 12 lanes, 2 (9534190031, 9534190041) have 6 lanes.

# Formatting columns. Saved them into the pheno_and_batch sheet in 
# Subject_Demographics_April2015_FINAL.ods

#######################################
# 1.2 the expression table 

mkdir input_files

# gs_report 'Sample_and_Control_Probe_Profile_FinalReport.txt'
cp Genome\ studio/BATCH3\ 91014/Groupprobes.txt \
    input_files/Sample_and_Control_Probe_Profile_FinalReport.txt
# basically a table of expression levels.
# 7 lines of batch header
# 2 lines of group header
# 4542  columns x 47231 rows 
# TargetID ProbeID
# then 
# AVG_Signal
# MIN_Signal
# MAX_Signal
# NARRAYS     # all 1
# ARRAY_STDEV # all NaN
# BEAD_STDERR
# Avg_NBEADS
# Detection Pval
# 8 values per subject
# Then 28 other columns including GI, SYMBOL, PROBE_ID etc
# (4542-28-2) / 8 = 564 subjects

# some samples have 0 expression values.
# 9464921139_L 9402444035_B 9402444030_B 9402444017_B 
# correlation coefficents of them are 'NA'
# One sample 9402444009_B has negative correlation coefficents
# It has mostly 0 expression values.
# 
# 9402444009_B is the 122nd sample
head -n 20 input_files/Sample_and_Control_Probe_Profile_FinalReport.txt \
 | tail -n +9 | \
  awk -F"\t" '{i=122*8-6; print $(i+1),$(i+2),$(i+3),$(i+4),\
                                $(i+5),$(i+6),$(i+7),$(i+8)}'

# Put the 5 ids into the file to_remove_sample.ls.

# 9464921139_L
# 9402444035_B
# 9402444030_B
# 9402444017_B
# 9402444009_B

remove_sample.py > t.out 

mv t.out input_files/Sample_and_Control_Probe_Profile_FinalReport.txt

# gs_sample "sample_table_Final_Report.txt"
cp Genome\ studio/BATCH3\ 91014/sampletable.txt \
    input_files/sample_table_Final_Report.txt

#######################################
# 1.3 the batch file 

# tech_pheno_file 'batch_info.txt'
# get info from Subject_Demographics_April2015_FINAL.ods

# replace 'Unknown' and 'unknown' with 'NA'
# to change the format of dates
sed -i 's/\ January\ /\/01\//g' t.in
sed -i 's/\ February\ /\/02\//g' t.in
sed -i 's/\ March\ /\/03\//g' t.in
sed -i 's/\ April\ /\/04\//g' t.in
sed -i 's/\ May\ /\/05\//g' t.in
sed -i 's/\ June\ /\/06\//g' t.in
sed -i 's/\ July\ /\/07\//g' t.in
sed -i 's/\ August\ /\/08\//g' t.in
sed -i 's/\ September\ /\/09\//g' t.in
sed -i 's/\ October\ /\/10\//g' t.in
sed -i 's/\ November\ /\/11\//g' t.in
sed -i 's/\ December\ /\/12\//g' t.in
# use only the 559 ids

#######################################
# 1.4 the pheno file
# pheno_file 'pheno_info.txt'
# get info from Subject_Demographics_April2015_FINAL.ods

# "Sample.ID","SEX","GROUPS","TISSUE","PHENOTYPE","Study_ID"
# use only the 559 ids

#######################################
# 1.5 the controlprobe file

get_control.py > t.out
cp t.out input_files/controlprobe.txt

################################################################################
# 2. Microarry pre-processing workflow for Illumina BeadArray data             #
################################################################################

# https://github.com/KHP-Informatics/chip_gx/blob/master/\
# Illumina_expression_workflow/\
# GAP_illumina_gene_expression_workflow_preProcessing.md

# install  url zip xml2 etc devel packages
# in sudo R 
# source("http://bioconductor.org/biocLite.R")
# biocLite() the libraries

# library flashClust needed to be installed.

cd ba3 

# ba3.R in R
# source(paste(path_to_scripts, "/sjnewhouse_misc_R.R", sep = ""))

################################################################################
# 3. CSF measures to AD case/control classification                            #
################################################################################

# input: Subject_Demographics_with_chip_data_for_processing_April2015_FINAL.xml
# output: a classifier to tell AD or control from CSF measures

################################################################################
# 3.1 target

# column S (Status) of the table:
# 285 AD
#  64 CTL
#  10 Fronto-temporal lobe dementia
#   1 HD
#  11 HD_Early
#   7 HD_Moderate
# 141 MCI
#   2 Possible dementia with Lewy body
#   3 Possible vascular dementia
#   2 Probable dementia with Lewy body
#   4 Probable vascular dementia

# leave AD as AD, CTL as CTL, everything else as MCI
# 285 AD, 64 CTL, 181 MCI

# table elements AY86:AZ207 should be AY355:AZ476
################################################################################
# 3.2 input

# we have 205 CSF measures from the table 
# After removal of the duplications, left with 113 samples in csf_in.ods
# column AW, AY and BA

#             ab42              pTau             tTau
# range:      65.250~776.045    10.86~104.105    35.775~367.345  
# mean:       386               37.9             115 
# sdv:        141               19.4             66.5
# normal (0): 52 (<390)         62 (<35)         52 (<98) 
# case (1):   61 (>390)         51 (>35)         61 (>98)  

#       2 AD  0 0 0     6 CTL  0 0 0
#       
#       3 AD  1 0 0     4 CTL  1 0 0
#                       2 CTL  0 1 0
#       2 AD  0 0 1           
#         
#       5 AD  0 1 1           
#       8 AD  1 0 1           
#       
#      16 AD  1 1 1           
#
#      19 MCI 0 0 0
#       6 MCI 0 0 1
#       4 MCI 0 1 0
#       6 MCI 0 1 1
#      10 MCI 1 0 0
#       2 MCI 1 0 1
#       2 MCI 1 1 0
#      16 MCI 1 1 1
#

# After bpnn.py single hidden unit back-propagation network training,
# the formula is 
#    ab42/188.7 - pTau/43 - tTau/58 + 0.0921 > 0 (CTL)
# or ab42/188.7 - pTau/43 - tTau/58 + 0.0921 < 0 (AD)

################################################################################
# validation of ba2 random forest                                              #
################################################################################

# 1. probes 
http://www.j-alz.com/issues/33/vol33-3.html#supplementarydata03

# 2. software 
http://www.imbs-luebeck.de/imbs/taxonomy/term/1
# It asks for Boost and gsl. 
# switching function locations in TermResult.h to avoid error in compiling.

# the training script from Martina
# rjungle --file=train_set.dat --treetype=1 --ntree=500 --mtry=500 \
#     -B 3 --impmeasure=5 --nimpvar=100 --memmode=0 --depvarname=Class \
#     --seeed=2178 --outprefix=train.result 

# What I did with 200 randomly chosen SNPs from pe19 set

plink --noweb --bfile pe19 --extract t2.ls --recodeA
cut -d' ' -f6- plink.raw  > t.out
sed -i 's/-9/NA/' t.out 
mv t.out plink_raw.data

rjungle --file=plink_raw.data \
    --treetype=1 --ntree=100 \
    --mtry=100 -B 3 --impmeasure=5 --nimpvar=100 --memmode=0 \
    --depvarname=PHENOTYPE --seeed=2178 --outprefix=train.result \
    --impute=30
# takes a few hours
# didn't give me the trained trees? 

