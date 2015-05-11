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


 