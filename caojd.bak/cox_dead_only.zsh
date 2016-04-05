################################################################################
# COX analysis on Dead only                                                    #
################################################################################

cd /home/ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19
mkdir GenABEL_DeadOnly 

################################################################################
# 1. convert the files                                                         #
################################################################################

# To run Cox analysis, we need 4 types of files:
# 1. phenotype file 
# 2. dosage file   (in MACH dosage format or DatABEL fvi format)
# 3. SNP info file (got in the previous COX analysis)
# 4. SNP map file  (got in the previous COX analysis)

for cohort in D1 D2 IR MGH NIT SLA UK ; do
    awk '{if ($NF == "1") print $1}' ONLY_DEAD_REG/${cohort}_pheno.sample \
        > ${cohort}_do.ls 
done 

rm t.out 
for cohort in D1 D2 IR MGH NIT SLA UK ; do
    grab -f ${cohort}_do.ls  ONLY_DEAD_REG/${cohort}_DEAD.sample | \
        awk '{print $NF}' >> t.out 
done 




tail -n +2 ONLY_DEAD_REG/IR_pheno.sample | \
    awk '{if ($NF==1) print $1}' > IR_do.ls 



#######################################
# 1.1 phenotype file 

# sample files:
for cohort in D1 D2 IR MGH NIT SLA UK ; do
    wc ONLY_DEAD_REG/${cohort}_pheno.sample
done 

# check ids have no duplications 
for cohort in D1 D2 IR MGH NIT SLA UK ; do
    tail -n +2 ONLY_DEAD_REG/${cohort}_pheno.sample \
        | awk '{print $1}' | sort | uniq -d 
done 

# check ids are in the complete pheno files:
for cohort in D1 D2 IR MGH NIT SLA UK ; do
    tail -n +2 ONLY_DEAD_REG/${cohort}_pheno.sample | \
        awk '{print $1}' > t.ls 
    grab -f t.ls GenABEL_OUTPUT/${cohort}_1.pheno > t.out 
    wc t.ls t.out 
done 
# got problem with SLA again. '1','2','4' will not match 'A1','A2','A115'


# get the lists of ids of dead only 
for cohort in D1 D2 IR MGH NIT SLA UK ; do
    tail -n +2 ONLY_DEAD_REG/${cohort}_pheno.sample \
        | awk '{print $1}' > ${cohort}_do.ls
done 

tail -n +2 ONLY_DEAD_REG/SLA_pheno.sample | awk '{print $2}' > SLA_do.ls 
sort SLA_do.ls | uniq -d 

# four types of analysis:
# different covariates
# 1. pc, sex and age 
# 2. pc and sex 
# 3. pc and age
# 4. pc 


for cohort in D1 D2 IR MGH NIT SLA UK ; do
    head -n 1 GenABEL_INPUT/${cohort}_pheno  > t.head 
    grab -f ${cohort}_do.ls GenABEL_INPUT/${cohort}_pheno > t.out 
    cat t.head t.out > GenABEL_DeadOnly/${cohort}_1.pheno

    awk        '{$5="";print $0}' GenABEL_DeadOnly/${cohort}_1.pheno \
        > GenABEL_DeadOnly/${cohort}_2.pheno
    
    awk        '{$4="";print $0}' GenABEL_DeadOnly/${cohort}_1.pheno \
        > GenABEL_DeadOnly/${cohort}_3.pheno
    
    awk '{$5="";$4=""; print $0}' GenABEL_DeadOnly/${cohort}_1.pheno  \
        > GenABEL_DeadOnly/${cohort}_4.pheno
done

for cohort in D1 D2 IR MGH NIT SLA UK ; do 
    awk '{print $3}' GenABEL_DeadOnly/${cohort}_1.pheno | sort | uniq -c ;
done

#######################################
# 1.2 dosage files

# 1.2.1 gen files

# after got all the bgen files 

cd /home/ifoghbrc/scratch/BUILD_37_ALS_GWAS/COX_h19/GenABEL_DeadOnly


get_mach.R
#--------------------------------------
library(MASS,lib='~/scratch/temp/R-2.15.2/library')
library(GenABEL,lib='~/scratch/temp/R-2.15.2/library')
library(DatABEL,lib='~/scratch/temp/R-2.15.2/library')

args = commandArgs(trailingOnly = TRUE)

cohort=args[1]
chr=args[2]

in_gfile_name=paste(cohort,'_',chr,'.gen', sep='')
in_ifile_name=paste(cohort,'_',chr,'.info',sep='')
in_sfile_name=paste(cohort,'_dead_only.sample',sep='')

out_mach_name=paste('mach_',cohort,'_',chr,sep='')

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

#######################################
# get gen file 

/share/bin/qctool -omit-chromosome \
    -g   ${cohort}_deal_only_${chr}.bgen \
    -og  $t_{cohort}_${chr}.gen

grep -v MERGED_DEL $t_{cohort}_${chr}.gen > ${cohort}_${chr}.gen 

rm  $t_{cohort}_${chr}.gen

#######################################
# convert into mach dosage files

~/scratch/temp/R-2.15.2/bin/Rscript get_mach.R ${cohort} ${chr}

#######################################
# remove temp files 

rm ${cohort}_${chr}.gen 
rm ${cohort}_${chr}.gen.prob.* 

rm mach_${cohort}_${chr}.machdose

#--------------------------------------
# the end of cox_convert.sh 


for cohort in D1 D2 IR MGH NIT SLA UK ; do
    for chr in 01 02 03 04 05 06 07 08 09 10 11 \
               12 13 14 15 16 17 18 19 20 21 22 ; do
        qsub cox_convert.sh $cohort $chr 
    done 
done

