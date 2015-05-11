################################################################################
#                             Cox analysis                                     #
################################################################################

################################################################################
# 1. install the software packages                                             #
################################################################################

#######################################
# 1.1 ProbABEL

# ProbABEL version 0.4.1 was installed last time (Sept 2013) at 
~/scratch/temp/probabel-0.4.1/

# It is working properly. No need to install a newer version. 
# The newer version (0.4.3) will not work on the cluster because of the 
# lacking of libeigen3 library. 
# However, the older version (0.4.1) is fine with that.

# In case you need to install it again:
cd ~/scratch/temp
wget http://www.genabel.org/sites/default/files/software/probabel-0.4.1.tar.gz
tar xvzf probabel-0.4.1.tar.gz
cd probabel-0.4.1
./configure --without-eigen 
make 
# need to press the 'Enter' key a few times when 'Latex' runs. 
cd ..

# now the program 'pacoxph' is available at
~/scratch/temp/probabel-0.4.1/src/pacoxph

#######################################
# 1.2 GenABEL 

# We need the MACH dosage files and DatABEL filevector files  
# for the Cox analysis. 

# We need the GenABEL R package. 
# The combination of R 3.0.2 and GenABEL 1.8-0 on the cluster doesn't work. 
# The author promised to fix it about one year ago. Won't wait for that. 

# We need an older R (2.15.2) and an older GenABEL (1.7-3).

###################
# 1.2.1 install R 2.15.2

cd ~/scratch/temp

wget http://cran.cnr.berkeley.edu/src/base/R-2/R-2.15.2.tar.gz

tar xvzf R-2.15.2.tar.gz 

cd R-2.15.2/
./configure
make
# It takes quite a while to make. 
cd ..

# now the older R is available at 
~/scratch/temp/R-2.15.2/bin/R 

###################
# 1.2.2 install GenABEL 1.7-3

wget ftp://cran.r-project.org/pub/R/src/contrib/\
Archive/GenABEL/GenABEL_1.7-3.tar.gz
wget ftp://cran.r-project.org/pub/R/src/contrib/\
Archive/DatABEL/DatABEL_0.9-2.tar.gz

~/scratch/temp/R-2.15.2/bin/R 
# in R

install.packages('GenABEL_1.7-3.tar.gz', repos = NULL, type="source",
	             lib='~/scratch/temp/R-2.15.2/library')

install.packages('DatABEL_0.9-2.tar.gz',repos = NULL, type="source",
	             lib='~/scratch/temp/R-2.15.2/library')

# now GenABEL can be loaded with 
~/scratch/temp/R-2.15.2/bin/R 

library(MASS,lib='~/scratch/temp/R-2.15.2/library')
library(GenABEL,lib='~/scratch/temp/R-2.15.2/library')

#######################################
# 1.3 cleaning
cd ~/scratch/temp
rm probabel-0.4.1.tar.gz
rm R-2.15.2.tar.gz
rm GenABEL_1.7-3.tar.gz
rm DatABEL_0.9-2.tar.gz

################################################################################
# 2. convert the files                                                         #
################################################################################

# To run Cox analysis, we need 4 types of files:
# 1. phenotype file 
# 2. dosage file   (in MACH dosage format or DatABEL fvi format)
# 3. SNP info file (in MACH info format)
# 4. SNP map file  (in MACH legend format)

cd ~/scratch/BUILD_37_ALS_GWAS/COX_h19/

# info, bgen and pheno files in 
GenABEL_INPUT/
# sample files in
./
# output files to 
GenABEL_OUTPUT/


get_mach.R
#--------------------------------------
library(MASS,lib='~/scratch/temp/R-2.15.2/library')
library(GenABEL,lib='~/scratch/temp/R-2.15.2/library')
library(DatABEL,lib='~/scratch/temp/R-2.15.2/library')

args = commandArgs(trailingOnly = TRUE)

out_dir=args[1]
cohort=args[2]
chr=args[3]

in_gfile_name=paste(out_dir,'/',cohort,'_',chr,'.gen', sep='')
in_ifile_name=paste(out_dir,'/',cohort,'_',chr,'.info',sep='')
in_sfile_name=paste(cohort,'_CS.sample',sep='')

out_mach_name=paste(out_dir,'/mach_',cohort,'_',chr,sep='')

impute2mach(in_gfile_name,in_ifile_name,
            in_sfile_name,out_mach_name,maketextdosefile=TRUE)

#--------------------------------------
# the end of get_mach.R

cox_convert.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N cox_convert
#$ -cwd
#$ -q long.q,bignode.q,short.q
##$ -l h_vmem=5G

cohort=$1
chr=$2
# 01 02 ... 22

in_dir=GenABEL_INPUT
out_dir=GenABEL_OUTPUT

#######################################
# 1. get the gen file 
/share/bin/qctool -omit-chromosome \
    -g   ${in_dir}/${cohort}_CS-NO-DUPL_${chr}.bgen \
    -og ${out_dir}/t_${cohort}_${chr}.gen

#######################################
# 2. remove the 'MERGED_DEL' lines, otherwise getting Segmentation fault 
# with pacoxph
grep -v MERGED_DEL ${out_dir}/t_${cohort}_${chr}.gen \
                 > ${out_dir}/${cohort}_${chr}.gen 

rm ${out_dir}/t_${cohort}_${chr}.gen    

grep -v MERGED_DEL ${in_dir}/${cohort}_${chr}.no-dup.info \
                  > ${out_dir}/${cohort}_${chr}.info 

#######################################
# 3. format into dosage files 
# need .dose.fvi, .machinfo and .machlegend files
# all are relatively small 

~/scratch/temp/R-2.15.2/bin/Rscript get_mach.R ${out_dir} ${cohort} ${chr}

#######################################
# 4. remove the temp files 

rm ${out_dir}/${cohort}_${chr}.gen 
rm ${out_dir}/${cohort}_${chr}.gen.prob.* 

rm ${out_dir}/mach_${cohort}_${chr}.machdose

#--------------------------------------
# the end of cox_convert.sh

for cohort in D1 D2 IR MGH NIT SLA UK ; do
    for chr in 01 02 03 04 05 06 07 08 09 10 11 \
               12 13 14 15 16 17 18 19 20 21 22 ; do
        qsub cox_convert.sh $cohort $chr 
    done 
done

################################################################################
# 3. run Cox analysis                                                          #
################################################################################

# four types of analysis:
# different covariates
# 1. pc, sex and age 
# 2. pc and sex 
# 3. pc and age
# 4. pc 

cd ~/scratch/BUILD_37_ALS_GWAS/COX_h19/

#######################################
# 3.1 prepare pheno tables

for cohort in D1 D2 IR MGH NIT SLA UK ; do
    cp GenABEL_INPUT/${cohort}_pheno GenABEL_OUTPUT/${cohort}_1.pheno 
    
    awk        '{$5="";print $0}' GenABEL_INPUT/${cohort}_pheno \
        > GenABEL_OUTPUT/${cohort}_2.pheno
    
    awk        '{$4="";print $0}' GenABEL_INPUT/${cohort}_pheno \
        > GenABEL_OUTPUT/${cohort}_3.pheno
    
    awk '{$5="";$4=""; print $0}' GenABEL_INPUT/${cohort}_pheno \
        > GenABEL_OUTPUT/${cohort}_4.pheno
done
 
#######################################
# 3.2 run the analysis

cox_run.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N cox_run
#$ -cwd
#$ -q long.q,bignode.q,short.q
##$ -l h_vmem=5G

out_dir=GenABEL_OUTPUT

cohort=$1
chr=$2
# 01 02 ... 22
chr_int=$(echo $chr | sed 's/^[0]//')
# 1 2 ... 22

~/scratch/temp/probabel-0.4.1/src/pacoxph \
    -p ${out_dir}/${cohort}_1.pheno   \
    -d ${out_dir}/${cohort}_${chr}.gen.dose.fvi \
    -i ${out_dir}/mach_${cohort}_${chr}.machinfo \
    -m ${out_dir}/mach_${cohort}_${chr}.machlegend \
    -c ${chr_int} \
    -o ${out_dir}/${cohort}_${chr}_1_cox 

~/scratch/temp/probabel-0.4.1/src/pacoxph \
    -p ${out_dir}/${cohort}_2.pheno   \
    -d ${out_dir}/${cohort}_${chr}.gen.dose.fvi \
    -i ${out_dir}/mach_${cohort}_${chr}.machinfo \
    -m ${out_dir}/mach_${cohort}_${chr}.machlegend \
    -c ${chr_int} \
    -o ${out_dir}/${cohort}_${chr}_2_cox 

~/scratch/temp/probabel-0.4.1/src/pacoxph \
    -p ${out_dir}/${cohort}_3.pheno   \
    -d ${out_dir}/${cohort}_${chr}.gen.dose.fvi \
    -i ${out_dir}/mach_${cohort}_${chr}.machinfo \
    -m ${out_dir}/mach_${cohort}_${chr}.machlegend \
    -c ${chr_int} \
    -o ${out_dir}/${cohort}_${chr}_3_cox 

~/scratch/temp/probabel-0.4.1/src/pacoxph \
    -p ${out_dir}/${cohort}_4.pheno   \
    -d ${out_dir}/${cohort}_${chr}.gen.dose.fvi \
    -i ${out_dir}/mach_${cohort}_${chr}.machinfo \
    -m ${out_dir}/mach_${cohort}_${chr}.machlegend \
    -c ${chr_int} \
    -o ${out_dir}/${cohort}_${chr}_4_cox 

#--------------------------------------
# the end of cox_run.sh 


for cohort in D1 D2 IR MGH NIT SLA UK ; do
    for chr in 01 02 03 04 05 06 07 08 09 10 11 \
               12 13 14 15 16 17 18 19 20 21 22 ; do
        qsub cox_run.sh $cohort $chr 
    done 
done


################################################################################
# 4 Cox analysis with site-onset and with deal-only                            #
################################################################################

cd /home/klinbrc/scratch/dev/cox

# in /home/ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_OUTPUT
# there are four types of COX output there: 

# 1. pc, sex and age 
# 2. pc and sex 
# 3. pc and age
# 4. pc 

# we need the fifth one: 5. pc, sex, age and site-onset 

# the site-onset info is available at 
~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/PHENO_INFO/${cohort}_site.sample

# the type 1 pheno files are at 
~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_OUTPUT/${cohort}_1.pheno

#######################################
# 4.1 to get type 5 pheno files

site_dir=~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/PHENO_INFO
in_dir=~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_OUTPUT

cohort=IR

awk '{print $1}' ${in_dir}/${cohort}_1.pheno | tail -n +2 > ind.ls 
# get the individue id list 

head -n 2 $site_dir/${cohort}_site.sample 

head -n 2 $site_dir/*_site.sample 
# the site-onset column are at different positions of the sample files 
# column 10 for most, column 7 for NIT

# ID_2 of the sample files should be used (expecially for SLA)
# IR has only NA 

t.head 
#--------------------------------------
site
#--------------------------------------
# end of t.head 


site_dir=~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/PHENO_INFO
in_dir=~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_OUTPUT

for cohort in D1 D2 IR MGH NIT SLA UK ; do
    site_col=10
    if [ "$cohort" == "NIT" ]; then 
        site_col=7
    fi
    # echo $site_col 
    
    in_sample_file=$site_dir/${cohort}_site.sample 
    
    awk -v sc=$site_col '{print $2,$sc}' $in_sample_file | \
        tail -n +2 > ${cohort}.sites 
    
    in_pheno_file=${in_dir}/${cohort}_1.pheno
    tail -n +2 $in_pheno_file | awk '{print $1}' > ind.ls 
    # get the individue id list 
    
    sort_table -f ind.ls ${cohort}.sites  > ${cohort}.out 
    awk '{print $2}'  ${cohort}.out  > t.out 
    cat t.head t.out > t.ls 
    # the site-onset list
    
    rm  ${cohort}.out  ${cohort}.sites 
    
    paste $in_pheno_file t.ls > ${cohort}_5.pheno 
done 

# check column numbers
for cohort in D1 D2 IR MGH NIT SLA UK ; do
    awk '{ print NF}'  ${cohort}_5.pheno  | sort -g | uniq -c 
done 

#######################################
# 4.2 run the cox analysis with site and all other covariates

cox_run.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N cox_run
#$ -cwd
#$ -q long.q,bignode.q,short.q
##$ -l h_vmem=5G

in_dir=~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_OUTPUT

cohort=$1
chr=$2
# 01 02 ... 22
chr_int=$(echo $chr | sed 's/^[0]//')
# 1 2 ... 22

~/scratch/temp/probabel-0.4.1/src/pacoxph \
    -p ${cohort}_5.pheno   \
    -d ${in_dir}/${cohort}_${chr}.gen.dose.fvi \
    -i ${in_dir}/mach_${cohort}_${chr}.machinfo \
    -m ${in_dir}/mach_${cohort}_${chr}.machlegend \
    -c ${chr_int} \
    -o ${cohort}_${chr}_5_cox 
#--------------------------------------
# end of cox_run.sh 

# It will NOT work on IR because of the all 'NA' site column.
# It will say: mematrix(): number of rows smaller then 1

# IR type 5 output should be identical to type 1 
# copy them 

in_dir=~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_OUTPUT

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    cp ${in_dir}/IR_${chr}_1_cox_add.out.txt IR_${chr}_5_cox_add.out.txt
done 

# run cox analysis for the other 6 cohorts

for cohort in D1 D2 MGH NIT SLA UK ; do
    for chr in 01 02 03 04 05 06 07 08 09 10 11 \
               12 13 14 15 16 17 18 19 20 21 22 ; do
        qsub cox_run.sh ${cohort} ${chr}
    done 
done 

mkdir full_set_type5
mv *_5_cox_add.out.txt full_set_type5/

#######################################
# 4.3 Dead only sets

cd /home/klinbrc/scratch/dev/cox

# get the pheno files 
cp ~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_OUTPUT/*pheno . 

in_dir=~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_DeadOnly

for cohort in D1 D2 IR MGH NIT SLA UK ; do
	tail -n +3 $in_dir/${cohort}_dead_only.sample | awk '{print $2}' > ind.ls 
	# get the individual list 
	
	for type in 1 2 3 4 5 ; do 
	    head -n 1 ${cohort}_${type}.pheno  > t.head 
	    grab -f ind.ls ${cohort}_${type}.pheno > t.out 
	    cat t.head t.out > ${cohort}_${type}.pheno
	done 
done 
# finish getting pheno files 
# Again, IR type 5 pheno file has an all 'NA' site column. 
# SLA has 150 'NA' out of 201 individuals, not good either. 

cox_run.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N cox_run
#$ -cwd
#$ -q long.q,bignode.q,short.q
##$ -l h_vmem=5G

in_dir=~ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_DeadOnly

cohort=$1
chr=$2                                     # 01 02 ... 22
chr_int=$(echo $chr | sed 's/^[0]//')      #  1  2 ... 22

for type in 1 ; do 
    ~/scratch/temp/probabel-0.4.1/src/pacoxph \
        -p  ${cohort}_${type}.pheno   \
        -d ${in_dir}/${cohort}_${chr}.gen.dose.fvi \
        -i ${in_dir}/mach_${cohort}_${chr}.machinfo \
        -m ${in_dir}/mach_${cohort}_${chr}.machlegend \
        -c ${chr_int} \
        -o  ${cohort}_${chr}_${type}_cox 
done 

#--------------------------------------
# end of cox_run.sh 

for cohort in D1 D2 MGH NIT SLA UK ; do
    for chr in 01 02 03 04 05 06 07 08 09 10 11 \
               12 13 14 15 16 17 18 19 20 21 22 ; do
        qsub cox_run.sh ${cohort} ${chr}
    done 
done 

# run IR differently, without type 5 

# 13 sets of genotype missing
# rebuild the genotype files 

# the error messages are some rsid are not unique 

qsub cox_convert.sh D1  15
qsub cox_convert.sh D2  06
qsub cox_convert.sh D2  15
qsub cox_convert.sh MGH 15
qsub cox_convert.sh NIT 01
qsub cox_convert.sh NIT 15
qsub cox_convert.sh SLA 01
qsub cox_convert.sh SLA 06
qsub cox_convert.sh SLA 15
qsub cox_convert.sh UK  01
qsub cox_convert.sh UK  15
qsub cox_convert.sh IR  01
qsub cox_convert.sh IR  15    rs1000290 



mkdir dead_only 
mv *add.out.txt dead_only 

################################################################################
# the end of Cox analysis                                                      #
################################################################################
