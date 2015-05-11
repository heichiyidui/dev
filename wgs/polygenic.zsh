################################################################################
# polygenic score calculation                                                  #
################################################################################

################################################################################
# input:                                                                       #
# 1. SNP score table: 130K SNPs, Allele, OR                                    #
# 2. GEN files                                                                 #
# 3. sample files                                                              #
################################################################################

# SNP score table
# GERAD_pruned025_summary_stats.txt

# sample files:
# /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/SAMPLE_files

# gen files:
# /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/*/*Imputed_Chromosomes/

################################################################################
# steps:                                                                       #
# 1. select subsets from GEN files                                             #
# 2. combine cohorts and tranform into plink style files                       #
# 3. calculate the polygenic scores                                            #
################################################################################

################################################################################
# 1. select subsets from GEN files                                             #
################################################################################

# working directory:
cd /home/klinbrc/scratch/temp

################################################################################
# 1.1 to get 22 snp lists

for chr in 22; do
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    chr_int=$(echo $chr | sed 's/^[0]//') # from 01 02 03 ... to 1 2 3 ...
    awk -v var="$chr_int" '{if ($2==var) print $1}' \
        GERAD_pruned025_summary_stats.txt | sort | uniq > snp_${chr}.ls
done

################################################################################
# 1.2 to get 22 gen files for each cohort

cohort="COR"
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    grab -f snp_${chr}.ls -c 2  /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/\
CORIELL_IMPUTATION/COR_imputed_chromosomes/\
NIH_Corriel_column_chr_${chr}.imputed.gen \
        > ${cohort}_${chr}.gen
done

cohort="D1"
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    grab -f snp_${chr}.ls -c 2 /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/\
DUTCH_1_IMPUTATION/D1_Impute_Chromosomes/Dutch1_column_chr_${chr}.imputed.gen \
        > ${cohort}_${chr}.gen
done

cohort="D2"
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    grab -f snp_${chr}.ls -c 2 /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/\
DUTCH_2_IMPUTATION/D2_Imputed_Chromosomes/Dutch2_column_chr_${chr}.imputed.gen \
        > ${cohort}_${chr}.gen
done

cohort="IR"
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    grab -f snp_${chr}.ls -c 2 /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/\
IRISH_IMPUTATION/IR_imputed_chromosomes/IRISH_column_chr_${chr}.imputed.gen \
        > ${cohort}_${chr}.gen
done

cohort="UK"
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    grab -f snp_${chr}.ls -c 2 /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/\
MERGE_UK_IMPUTATION/MUK_Chromosomes/UK_MERGE_chr_${chr}.gen \
        > ${cohort}_${chr}.gen
done

cohort="MGH"
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    grab -f snp_${chr}.ls -c 2 /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/\
META3_IMPUTATION/MGH_imputed_chromosomes/META3_column_chr_${chr}.imputed.gen \
        > ${cohort}_${chr}.gen
done

cohort="NIT"
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    grab -f snp_${chr}.ls -c 2 /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/\
NIT_NEW_IMPUTATION/NIT_Imputed_Chromosomes/NIH_IT_chr_${chr}.gen \
        > ${cohort}_${chr}.gen
done

cohort="SLA"
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    grab -f snp_${chr}.ls -c 2 /scratch/home/ifoghbrc/INTERNAT_IMPUTATION/\
SLAGEN_IMPUTATION/SLA_Imputed_Chromosomes/SLAGEN_column_chr_${chr}.imputed.gen \
        > ${cohort}_${chr}.gen
done

# around 5 mins or less for one chr gen file

#######################################
# A/T and C/G SNPS

tail -n +2 GERAD_pruned025_summary_stats.txt | awk '{print $4 $5}' | \
    sort | uniq -c

#  6138 AC
# 28004 AG
#    93 AT
#  5782 CA
#   140 CG
# 25343 CT
# 24721 GA
#   154 GC
#  5688 GT
#    83 TA
# 27861 TC
#  6171 TG

# only 93 + 140 + 154 + 83 A/T and C/G SNPs
# 0.36% of 130178 SNPs in the table
# can be removed

awk '{aa=$4 $5; if (aa=="AT" || aa=="TA" || aa=="CG" || aa=="GC") print $1}' \
  GERAD_pruned025_summary_stats.txt > to_remove_snp.ls
# 470 SNPs to remove

# for each cohort, combine chrs into one gen file
# we also need to removed duplicated SNPs

t1.py
#---------------------------------------
#!/home/klinbrc/bin/python3
import sys

to_remove_snps=open('to_remove_snp.ls').read().split()
to_remove_snps=set(to_remove_snps)

cohort=sys.argv[1] # 'COR', 'UK', 'IR', etc...
out_file=open(cohort+'.gen','w')
written_snps=set()

for chr in ['01','02','03','04','05','06','07','08','09','10','11',\
            '12','13','14','15','16','17','18','19','20','21','22']:
    in_file=open(cohort+'_'+chr+'.gen')
    for line in in_file:
        snp_id=line[:50].split()[1]
        if snp_id in to_remove_snps:
            continue
        if snp_id in written_snps:
            continue
        out_file.write(line)
        written_snps.add(snp_id)
    in_file.close()
out_file.close()
#---------------------------------------

for cohort in COR D1 D2 IR UK MGH NIT SLA ; do
    t1.py $cohort
done

# now the chr gen files can be removed
rm *_??.gen
rm snp_??.ls
rm to_remove_snp.ls
# we end up with 8 gen files for 8 cohorts

################################################################################
# 2. combine cohorts and tranform into plink style files                       #
################################################################################

################################################################################
# 2.1 transform gen and sample files into plink files

sdir=/scratch/home/ifoghbrc/INTERNAT_IMPUTATION/SAMPLE_files

/share/bin/gtool -G      --g COR.gen --s $sdir/CORIELL_AAO.sample \
    --phenotype STATUS --ped COR.ped     --map COR.map 

/share/bin/gtool -G      --g D1.gen  --s $sdir/DUTCH1_AAO.sample \
    --phenotype STATUS --ped D1.ped      --map D1.map 

/share/bin/gtool -G      --g D2.gen  --s $sdir/DUTCH2_AAO.sample \
    --phenotype STATUS --ped D2.ped      --map D2.map 

/share/bin/gtool -G      --g IR.gen  --s $sdir/IR_AAO.sample \
    --phenotype STATUS --ped IR.ped      --map IR.map 

/share/bin/gtool -G      --g UK.gen  --s $sdir/MERGE_UK_AAO.sample \
    --phenotype STATUS --ped UK.ped      --map UK.map 

/share/bin/gtool -G      --g MGH.gen --s $sdir/MGH_AAO.sample \
    --phenotype STATUS --ped MGH.ped     --map MGH.map 

/share/bin/gtool -G      --g NIT.gen --s $sdir/NIT_AAO.sample \
    --phenotype STATUS --ped NIT.ped     --map NIT.map 

/share/bin/gtool -G      --g SLA.gen --s $sdir/SLAGEN_AAO.sample \
    --phenotype STATUS --ped SLA.ped     --map SLA.map 

# not bothered about gender     

################################################################################
# 2.2 into plink binary files

for cohort in COR D1 D2 IR UK MGH NIT SLA ; do
	plink --noweb --file $cohort --make-bed --out $cohort --missing-genotype N 
done 

# fix the phenotypes 
for cohort in COR D1 D2 IR UK MGH NIT SLA ; do
   awk '{if ($6==1) $6=2; if ($6==-9) $6=1; print $0}' ${cohort}.fam > t.fam 
   mv t.fam ${cohort}.fam
done 

rm *.log *.nosex *.nof 
rm *.ped *.map 
rm *.gen 

################################################################################
# 2.3 fix SNP positions and strands

# update All SNP positions according to the SNP score table

t1.py
#--------------------------------------
#!/usr/bin/env python3
import sys

# read SNP position 
snp_position={}
score_table=open('GERAD_pruned025_summary_stats.txt')
score_table.readline()
for line in score_table:
    cols = line.split()
    snp_id  = cols[0]
    snp_position[snp_id] = cols[2]
score_table.close()

# update snp position in bim file 
ifile=open(sys.argv[1]) # input bim file 
for line in ifile:
    cols = line.split()
    snp_id = cols[1]
    cols[3] = snp_position[snp_id]
    print(' '.join(cols))
ifile.close()
#--------------------------------------

for cohort in COR D1 D2 IR UK MGH NIT SLA ; do
    t1.py ${cohort}.bim > t.bim 
    mv t.bim ${cohort}.bim
done 

# still found 8 C/G SNPs 
# 4590 with 0 in them such as 0C, 0G etc 
# remove them 

for cohort in COR D1 D2 IR UK MGH NIT SLA ; do
    rm to_rm_snp.ls 
    awk '{if ($5=="C" && $6=="G") print $2}' ${cohort}.bim >  to_rm_snp.ls 
    awk '{if ($5=="G" && $6=="C") print $2}' ${cohort}.bim >> to_rm_snp.ls 
    awk '{if ($5==0) print $2}'              ${cohort}.bim >> to_rm_snp.ls 

    plink --noweb --bfile $cohort --exclude to_rm_snp.ls \
        --make-bed --out $cohort
    rm to_rm_snp.ls 
done 

# find SNPs to flip according to the SNP score table

t1.py
#--------------------------------------
#!/usr/bin/env python3
import sys

# read SNP allele
snp_a1a2={}
score_table=open('GERAD_pruned025_summary_stats.txt')
score_table.readline() # the header line 
for line in score_table:
    cols = line.split()
    snp_id  = cols[0]
    snp_a1a2[snp_id] = cols[3] + cols[4]
score_table.close()

# update snp position in bim file 
ifile=open(sys.argv[1]) # input bim file 
for line in ifile:
    cols = line.split()
    snp_id = cols[1]
    a1 = cols[4]
    if a1 not in snp_a1a2[snp_id]:
        print(snp_id)
ifile.close()
#--------------------------------------

for cohort in COR D1 D2 IR UK MGH NIT SLA ; do
    t1.py ${cohort}.bim > to_flip_snp.ls 
    plink --noweb --bfile $cohort --flip to_flip_snp.ls --make-bed --out $cohort 
    rm to_flip_snp.ls 
done 

################################################################################
# 2.4 merge cohorts

plink --noweb --bfile COR --bmerge  D1.bed  D1.bim  D1.fam --make-bed --out m1
plink --noweb --bfile m1  --bmerge  D2.bed  D2.bim  D2.fam --make-bed --out m2
plink --noweb --bfile m2  --bmerge  UK.bed  UK.bim  UK.fam --make-bed --out m3
plink --noweb --bfile m3  --bmerge MGH.bed MGH.bim MGH.fam --make-bed --out m4
plink --noweb --bfile m4  --bmerge NIT.bed NIT.bim NIT.fam --make-bed --out m5
plink --noweb --bfile m5  --bmerge  IR.bed  IR.bim  IR.fam --make-bed --out m6

plink --noweb --bfile m6  --bmerge SLA.bed SLA.bim SLA.fam --make-bed --out m7
# got one SNP rs16922014 wrong in the last step
# A/G in m6, A/C in SLA
# remove it 

plink --noweb --bfile m6  --exclude m7.missnp --make-bed --out m6 
plink --noweb --bfile SLA --exclude m7.missnp --make-bed --out SLA 
plink --noweb --bfile m6  --bmerge SLA.bed SLA.bim SLA.fam --make-bed --out m7

# end up with a merged set of 13225 individuals, 129590 SNPs

##########################
# now clean the directory 
mv GERAD_pruned025_summary_stats.txt ..
mv m7.bim m7.fam m7.bed ..
rm * 
mv ../m7.??? .
mv ../GERAD_pruned025_summary_stats.txt .

################################################################################
# 3. calculate the polygenic scores                                            #
################################################################################

tail -n +2 GERAD_pruned025_summary_stats.txt | \
    awk '{print $1,$4,$6}' > ad_beta.raw 

plink --noweb --bfile m7 --score ad_beta.raw --out ad_beta
# plink thinks all individuals with missing gender are missing in phenotype too.
# assign all individuals with gender '1'
awk '{$5=1; print $0}' m7.fam > t.fam ; mv t.fam m7.fam 

plink --noweb --bfile m7 --score ad_beta.raw --out ad_beta
# of the 130178 SNP predictors in the raw score table, 129590 are found in m7

tar -czf m7.tgz m7.???
mv m7.tgz /home/klinbrc
mv ad_beta.profile /home/klinbrc

################################################################################
# 4. Regression                                                                #
################################################################################

R

poly=read.table('ad_beta.profile',header=TRUE)
fit <- lm(poly$PHENO ~ poly$SCORE)
summary(fit)

# Call:
# lm(formula = poly$PHENO ~ poly$SCORE)
#
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -0.4884 -0.4618 -0.4512  0.5379  0.5651 
#
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.462e+00  4.339e-03 336.867   <2e-16 ***
# poly$SCORE  1.276e+02  8.106e+01   1.574    0.115    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
# Residual standard error: 0.4985 on 13223 degrees of freedom
# Multiple R-squared:  0.0001874, Adjusted R-squared:  0.0001118 
# F-statistic: 2.479 on 1 and 13223 DF,  p-value: 0.1154

# p-value: 0.1154, not significant 

################################################################################
#                                                                              #
################################################################################
