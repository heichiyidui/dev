################################################################################
#                Finding variants in PCSK9 affecting LDL-C                     #
################################################################################

################################################################################
#                                                                              #
# The PCSK9 gene blocks LDL (Low-density lipoprotein) receptors. Less LDL      #
# receptor on the surface of liver to remove LDL from bloodstream, leads to    #
# higher LDL cholesterol (LDL-C) concentrations.                               #
# PCSK9 is a fat-control drug target. It is well studied with caucasian        #
# populations.                                                                 #
#                                                                              #
# Now verify our CKB dataset on the detection of SNPs associating with LDL-C   #
# level in this gene.                                                          #
#                                                                              #
################################################################################

# K:\kadoorie\Groups\Genetics\PROJECTS\PCSK9
# PCSK9\ analysis\ plan\ v3.docx

alias t='tail'
alias h='head'
alias skh='tail -n +2 '

################################################################################
# 1. genotype data                                                             #
################################################################################

################################################################################
# 1.1 the original set

# Start from the stage3 set at /kuser/shared/data/GWASphase12

plink --bfile /kuser/shared/data/GWASphase12/stage3 \
      --remove /kuser/shared/data/GWASphase12/stage3_mandatory_exclusions.txt \
      --autosome \
      --make-bed --out ckb_ph12_s3

# 636670 variants and 32205 people
# all with unknown father, mother and status
# no repeat in individual G-cryovial ids (ck_ids)

plink --bfile /kuser/shared/data/GWASphase12/stage3 \
      --remove /kuser/shared/data/GWASphase12/stage3_mandatory_exclusions.txt \
      --check-sex
# Obviously, the gender check has been done before.

################################################################################
# 1.2 very basic missing-call and MAF filter
# Do NOT perform it yet. Don't want miss anything yet.
# The plan says no QC on SNPs

# plink --bfile ckb_ph12_s3 \
#       --geno 0.05 \
#       --maf  0.0001 \
#       --make-bed --out ckb_ph12_s3_qc01

# # 51983 variants removed due to missing genotype data (--geno).
# # 30961 variants removed due to minor allele threshold(s)
# # 553726 variants and 32205 people left



################################################################################
# 1.3 heterogeneous related to missingness

# Do NOT remove any subjects yet.

# plink --bfile ckb_ph12_s3_qc01 --het
# # Not correct. plink says use a LD-free set.
# # A comparison was made.
# # The two sets of F values (from the whole and the LD-free sets)
# # are very close.
#
# plink --bfile ckb_ph12_s3_qc01 \
#       --geno 0.01 \
#       --hwe 1e-4 midp \
#       --maf 0.05 \
#       --indep-pairwise 1500 150 0.2
# # 120201 SNPs left
#
# plink --bfile ckb_ph12_s3_qc01 \
#       --extract plink.prune.in \
#       --make-bed --out t
#
# plink --bfile t  \
#       --het
#
# # and missingness
# plink --bfile ckb_ph12_s3_qc01 \
#       --missing

# paste plink.het plink.imiss | awk '{print $1,$2,$6,$12}' > t.in

# # in R
# #--------------------------------------
# library(ggplot2)
# data = read.table('t.in',header=T)
# p1 <- ggplot(data) + geom_point(aes(x=F, y=F_MISS)) + theme_bw()
# p2 <- ggplot(data=subset(data,F<0.10), aes(F)) +
#       geom_histogram(bins=30) + theme_bw()
#
# threshold = mean(data$F) - 3 * sd(data$F)
# data$lowF <- data$F < threshold
#
# # -0.0287 ~ 0.03048 for +-3 SD
# # 395 individuals flaged
# # or 20 individuals on the lower end
#
# p3 <- ggplot(data,aes(x=lowF, y=F_MISS)) + geom_boxplot( ) + theme_bw()
#
# write.table(subset(data,lowF)$IID,file='low_hom.ls',
#             row.names=F,col.names=F,quote=F)
# #--------------------------------------
# # end of the R script

# grab -f low_hom.ls -c 2 ckb_ph12_s3_qc01.fam | awk '{print $1}' | \
#     sort | uniq -c
# # The plate ids are scattered. I was expecting some plate contamination.

# # low homozygotes might mean sample contamination.
# # The 20 subjects with < 3 x SD F values are of very high missingness.
# # 0.022 vs 0.0029 for all subjects
# # And, according to the plink IBD estimation,
# # these subjects are related to almost EVERYONE else.
# # templed to remove them
# #######################################
# # don't remove any subject here yet

################################################################################
# 1.4 IBD and PCA

#######################################
# 1.4.1 LD pruning

plink --bfile ckb_ph12_s3 \
      --geno 0.01 \
      --maf 0.05 \
      --hwe 1e-4 midp \
      --indep-pairwise 1500 150 0.2

# 120201 SNPs in, 217496 out

plink --bfile ckb_ph12_s3 \
      --extract plink.prune.in \
      --make-bed --out pca

# 120201 variants and 32205 people

#######################################
# 1.4.2 IBD

plink --bfile pca \
      --genome --min 0.05

# It takes 20 minutes on NC2, uses many threads.
# plink.genome is huge without the min 0.05 filter

skh plink.genome | awk '{print $2 "\n" $4}' | sort | uniq > t.ls
# 32040 subjects to be removed? not acceptable.

skh plink.genome | awk '{print $2,$4}' > t.in

select_fam.py > to_remove.ls
# Remove (flag) the subjects according to their degrees of connection.
# The most related subjects are to be removed first.
# The top ones are still with low homozygotes.

# 6990 subjects flagged, 22% of 32205

# There are three pairs of people with PI_HAT about 0.5
# CK22775935 CK28153131
# CK28568808 CK28569055
# CK30580964 CK30589772
# The others are a huge family.

# Tried different ad hoc approaches to remove fewer nodes, but failed.
# Was doing the PCA with 7333 subjects flagged as related.
# Didn't change the results much.

# change the 6th column of pca.fam to 'rel' and 'no_rel'.

grab -f to_remove.ls -c 2 pca.fam | \
    awk '{print $1,$2,$3,$4,$5,"rel"}' > t.out
grab -f to_remove.ls -c 2 pca.fam -v | \
    awk '{print $1,$2,$3,$4,$5,"no_rel"}' >> t.out

awk '{print $2}' pca.fam > t.ls
sort_table -f t.ls -c 2 t.out > t2.out
mv t2.out pca.fam

printf "no_rel\n" > pop.ls

#######################################
# 1.4.3 PCA via EIGENSOFT

# By default, smartpca should be using multiple threads.
# First run PCA with the unrelated subjects only. Then project the factors to
# the related subjects.
# Don't remove outliers yet.
nohup /kuser/shared/bin/EIG/bin/smartpca.perl \
        -i pca.bed \
        -a pca.bim \
        -b pca.fam \
        -w pop.ls \
        -o no_rel.pca  \
        -p no_rel.plot \
        -e no_rel.eval \
        -l no_rel.log  \
        -m 0   &

# It took 13G of memory, about 24 hours on nc2.

# PCA on all subjects.
nohup /kuser/shared/bin/EIG/bin/smartpca.perl \
        -i pca.bed \
        -a pca.bim \
        -b pca.fam \
        -o all.pca  \
        -p all.plot \
        -e all.eval \
        -l all.log  \
        -m 0   &

# 21G of memory, about 48 hours on nc2


#######################################
# 1.4.4 PCA via plink (GCTA)

# PCA on unrelated subjects, then project back to the related.
awk '{print $1,$2,$6}' pca.fam > pca_cluster.in

nohup plink --bfile pca \
            --pca 10 \
            --within pca_cluster.in \
            --pca-cluster-names no_rel &

# plink uses about 80 threads
# 7G of memory
# about 5 hours

#######################################
# 1.4.5 check the PCA plots

# plot the eignsoft PCA results
ifile=no_rel.pca.evec
# or
# ifile=all.pca.evec

sed 's/\:/\t/' $ifile  -i

printf "CK24820387\nCK25228869\nCK28902540\nCK28730586\n" > t.ls
# found no study id for them. So there are no region codes either.

skh $ifile | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | \
    grep -f t.ls -v > t.in

# for ck_id.ls and study_id.ls, see section 2.
sort_table -f ck_id.ls t.in > t.out

printf "study_id rc pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 is_fam\n" > t.in

paste study_id.ls t.out | \
    awk '{print $1,substr($1,1,2),$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' \
        >> t.in

plot_pca.R

# plot the plink PCA results
awk '{print $6}' pca.fam > t.ls
paste plink.eigenvec t.ls > t.in

printf "CK24820387\nCK25228869\nCK28902540\nCK28730586\n" > t.ls
grep -v -f t.ls t.in > t.out

printf "study_id rc pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 is_fam\n" > t.in

paste study_id.ls t.out | \
    awk '{print $1,substr($1,1,2),$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,\
         $14}' >> t.in

plot_pca.R

# If we use the full set, after the first 2 PC2, the Suzhou family is draging
# other PCs everywhere.

# With relative-only PCA, we don't have such problem.

# PCs, up to PC7 and PC8, are obviously related to the region codes.
# plink is about 4 times faster than eigensoft.
# They produce mostly the same results.
#

################################################################################
# 1.5 Remove badly called SNPs found in the manual check.

# Do NOT do it yet.
# tail -n +2 ../ckb_batch_check/manual_chk_res.table | \
#     awk '{if ($2==0) print $1}' > to_remove_snp.ls
#
# # 13938 SNPs to be removed.

# # and remove the four subjects with missing ascertainments
#
# printf "CK24820387\nCK25228869\nCK28902540\nCK28730586\n" > t.ls
#
# grab -f t.ls -c 2 ckb_ph12_s3_qc02.fam > t.fam
#
# plink --bfile ckb_ph12_s3_qc02 \
#       --exclude to_remove_snp.ls \
#       --remove t.fam \
#       --make-bed --out ckb_ph12_s3_qc03
# # 546462 variants and 32181 people

################################################################################
# 2. phenotype data                                                            #
################################################################################

################################################################################
# 2.1 to get study ids

awk '{print $2}' pca.fam > ck_id.ls
# 32205 uniq ids in ck_id.ls

# 32410 subject ascerntaiments from
# GWAS_SNPdata_samples.xlsx in
# K:\kadoorie\Groups\Genetics\Data Archive\Project Sample Lists\Lists\

# The ids in the 'notes' column are absent from the ids obtained using the data
# request form. We removed this column.

# Some G-cryovial ids are of different format with the ids in the fam file.
# CK28185397-1 vs CK28185397-QC
# CK22754927-1 vs CK22754927-1-QC
# etc
# Changed the ids in the ascertainment file to confirm the plink fam files.

# 4 subjects in ckb_ph12_s3.fam are not found in the ascertainment files:
# CK24820387 CK25228869 CK28902540 CK28730586
# They were to be deleted from the genotype set?

# The modified and sorted file
GWAS_SNPdata_samples.csv

printf "CK24820387\nCK25228869\nCK28902540\nCK28730586\n" > t.ls
grep -f t.ls -v GWAS_SNPdata_samples.csv > t.out
mv t.out GWAS_SNPdata_samples.csv

skh  GWAS_SNPdata_samples.csv | awk -F"," '{print $1}' > ck_id.ls
skh  GWAS_SNPdata_samples.csv | awk -F"," '{print $2}' > study_id.ls

# 32201 uniq study and ck id pairs.
# ck_id.ls and study_id.ls were used in 1.5.1

################################################################################
# 2.2 stratification

# using the data request form to get snapshort 10 ldl data
# got 19423 numbers, scared by 1/100.

# 7114 got study_id of the genotype set
# The correlation coefficient is 0.35 very significant, but not high ...

#######################################
# the input files:
# ckb_ph12_s3.fam               # for ck_id in the genotype data set
# GWAS_SNPdata_samples.csv      # for corresponding study ids
# LDL-c_biochem_data.csv        # for direct LDL-C measures
# ind_ldl.csv                   # for indirect LDL-C measures
# age_base.csv                  # age at baseline
# age_resu1.csv                 # age at resurvey 1
# age_resu2.csv                 # age at resurvey 2
# PCSK9_sample_summary.csv      # sample ascertainment and survey id

get_strata.py > t.in
rint.R
# The table generated is in the file 't.out'.
# Maybe re-order the table according to the fam file?
awk '{print $2}' ckb_ph12_s3.fam > t.ls
sort_table -f t.ls t.out -c 2 -sh | grep -v NA > pheno.csv


# use plink no-rel 10 PCs and the pheno.csv for covariates
get_cov.py | sed 's/\ /\t/g' > cov.csv

################################################################################
# 3. linear regression and meta-analysis                                       #
################################################################################

################################################################################
# 3.1 linear regression

for stratum in {1..6}; do
    awk -vStr=$stratum '{if ($5==Str) print $2 }' pheno.csv > t.ls
    grab -f t.ls -c 2 ckb_ph12_s3.fam > st$stratum.fam
done

wc st?.fam
# 4762  28572 153184 st1.fam
# 5210  31260 167661 st2.fam
#  167   1002   5363 st3.fam
# 1265   7590  40848 st4.fam
# 6687  40122 215004 st5.fam
# 4177  25062 135136 st6.fam

#######################################
# 3.1.1 raw LDL-C measures

# for the big strata
for st in st1 st2 st4 st5 st6 ; do
    nohup plink --bfile ckb_ph12_s3 \
      --keep $st.fam \
      --pheno pheno.csv \
      --pheno-name ldl_c \
      --covar cov.csv \
      --linear hide-covar --ci 0.95 \
      --out $st.raw &
done

# for the smallest stratum, using only a few covariates.
# otherwise we will get all 'NA'

nohup plink --bfile ckb_ph12_s3 \
  --keep st3.fam \
  --pheno pheno.csv \
  --pheno-name ldl_c \
  --covar cov.csv \
  --covar-name sex,age,pc1,pc2 \
  --linear hide-covar --ci 0.95 \
  --out st3.raw &

#######################################
# 3.1.2 RINT LDL-C measures
# note the '--covar-name' lines are different

for st in st1 st2 st4 st5 st6 ; do
    nohup plink --bfile ckb_ph12_s3 \
      --keep $st.fam \
      --pheno pheno.csv \
      --pheno-name rint_ldl_c \
      --covar cov.csv \
      --covar-name pc1-pc10 \
      --linear hide-covar --ci 0.95 \
      --out $st.rint &
done

# for the smallest stratum, using only a few covariates.
nohup plink --bfile ckb_ph12_s3 \
  --keep st3.fam \
  --pheno pheno.csv \
  --pheno-name rint_ldl_c \
  --covar cov.csv \
  --covar-name pc1,pc2 \
  --linear hide-covar --ci 0.95 \
  --out st3.rint &

# RINT LDL-C can also be used in a cancatenated single set
nohup plink --bfile ckb_ph12_s3 \
    --pheno pheno.csv \
    --pheno-name rint_ldl_c \
    --covar cov.csv \
    --covar-name pc1-pc10 \
    --linear hide-covar --ci 0.95 \
    --out rint &

# Well, there's trouble.
# In the result, we found 4 APOE SNPs with P value 0.
# AX-59878593
# AX-32750083
# AX-82989693
# AX-94348016
# Changed all of them into 4.574e-292, the smallest P value otherwise.
# Or maybe the smallest double number greater than zero?
# 4.94066e-324?

################################################################################
# 3.2 QQ, Manhattan plots and Lambda

#######################################
# QQ and Manhattan plots

# let's have a look
for st in st1 st2 st3 st4 st5 st6 ; do
    plot_qq_man.R $st.raw.assoc.linear &
done

for st in st1 st2 st3 st4 st5 st6 ; do
    plot_qq_man.R $st.rint.assoc.linear &
done

plot_qq_man.R rint.assoc.linear &

#######################################
# to calculate lambda
# lambda.R
#--------------------------------------
#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
ifile_name = args[1]

data=read.table(ifile_name,header=T)
data=subset(data,!is.na(P))

chisq <- qchisq(1-data$P,1)
lambda = median(chisq)/qchisq(0.5,1)
cat(args[1],'lambda: ',lambda,'\n')
#--------------------------------------
# end of lambda.R

for st in st1 st2 st3 st4 st5 st6 ; do
    lambda.R $st.raw.assoc.linear
done

for st in st1 st2 st3 st4 st5 st6 ; do
    lambda.R $st.rint.assoc.linear
done

lambda.R rint.assoc.linear

#   stratum     size  lambda_raw  lambda_RINT
#   stratum 1	4762  1.007485    1.019263
#   stratum 2	5210  1.008424    1.010303
#   stratum 3	 167  0.9920914   0.9990672
#   stratum 4	1265  0.9907005   0.998135
#   stratum 5	6687  1.024951    1.026377
#   stratum 6	4177  1.019736    1.02971
#   all_rint   22268  NA          1.056657

################################################################################
# 3.3 meta-analysis

# put A2 into the tables
for ifile in st?.*.assoc.linear ; do
    add_a2.py $ifile > t.out
    mv t.out $ifile
done

# meta
metal < t_metal.sh

# try plot the meta results and calculate the Lambda
# first add the CHR SNP and BP columns

printf "CHR MarkerName BP\n" > t.loc
awk '{print $1,$2,$4}' ckb_ph12_s3.bim >> t.loc

ifile='ldl_direct1.tbl'
# ifile='ldl_all1.tbl'
awk '{print $1}' $ifile > t.ls
grab -f t.ls t.loc -c 2 > t.out
sort_table -f t.ls t.out -c 2 | sed 's/MarkerName/SNP/' >  t.in

paste t.in $ifile | sed 's/P-value/P/' > t.out
plot_qq_man.R t.out
lambda.R t.out

mv  t.out_qq.png  $ifile\_qq.png
mv  t.out_man.png $ifile\_man.png

# 1.040715 direct meta
# 1.050356 all meta

# again, we have 0 P values in all strata meta
# again, they are in the Apoe region
# AX-59878593  2.52e-334
# AX-32750083  2.09e-337
# The values are so small, they are considered 0 in R.
# manually change them to 2.52e-324

################################################################################
# 4. In the PCSK9 region, conditioned on leading SNPs                          #
################################################################################

# looking at the PCSK9 region only
plink --bfile ckb_ph12_s3 \
      --from AX-105169173 --to AX-31657601 \
      --make-bed --out geno

# 160 variates, 32205 subjects

################################################################################
# 4.1 RINT in a single set, plink only
# This is the simplest one?

printf "" > c_snp.ls
plink --bfile geno \
    --pheno pheno.csv \
    --pheno-name rint_ldl_c \
    --covar cov.csv \
    --covar-name pc1-pc10 \
    --linear hide-covar --ci 0.95 \
    --condition-list c_snp.ls \
    --out rint_p9_0

# the leading SNP is AX-83389438

skh rint_p9_0.assoc.linear | grep -v NA | sort -g -k 12 | \
    h -n 1 | awk '{print $2}' >> c_snp.ls

plink --bfile geno \
    --pheno pheno.csv \
    --pheno-name rint_ldl_c \
    --covar cov.csv \
    --covar-name pc1-pc10 \
    --linear hide-covar --ci 0.95 \
    --condition-list c_snp.ls \
    --out rint_p9_1

# the second SNP is AX-39912161
skh rint_p9_1.assoc.linear | grep -v NA | sort -g -k 12 | \
    h -n 1 | awk '{print $2}' >> c_snp.ls

plink --bfile geno \
    --pheno pheno.csv \
    --pheno-name rint_ldl_c \
    --covar cov.csv \
    --covar-name pc1-pc10 \
    --linear hide-covar --ci 0.95 \
    --condition-list c_snp.ls \
    --out rint_p9_2

# the third SNP is AX-11576926
skh rint_p9_2.assoc.linear | grep -v NA | sort -g -k 12 | \
    h -n 1 | awk '{print $2}' >> c_snp.ls

plink --bfile geno \
    --pheno pheno.csv \
    --pheno-name rint_ldl_c \
    --covar cov.csv \
    --covar-name pc1-pc10 \
    --linear hide-covar --ci 0.95 \
    --condition-list c_snp.ls \
    --out rint_p9_3
# No SNP after that
skh rint_p9_3.assoc.linear | grep -v NA | sort -g -k 12 | h

# The leading SNPs :
# AX-83389438
# AX-39912161
# AX-11576926

# If we removed AX-83389438, followed by AX-31642001,
# we would go on to remove AX-11541856
# and AX-31642169 after that.

for ifile in rint_p9_?.assoc.linear ; do
    plot_bp_p.R $ifile
done

################################################################################
# 4.2 raw plink, then raw metal, then all metal
printf "" > c_snp.ls

for st in st1 st2 st4 st5 st6; do
    plink --bfile geno \
      --keep $st.fam \
      --pheno pheno.csv \
      --pheno-name ldl_c \
      --covar cov.csv \
      --linear hide-covar --ci 0.95 \
      --condition-list c_snp.ls \
      --out $st.raw
done

plink --bfile geno \
  --keep st3.fam   \
  --pheno pheno.csv \
  --pheno-name ldl_c  \
  --covar cov.csv  \
  --covar-name sex,age,pc1,pc2 \
  --linear hide-covar --ci 0.95 \
  --condition-list c_snp.ls \
  --out st3.raw

for st in st1 st2 st3 st4 st5 st6 ; do
    add_a2.py $st.raw.assoc.linear > t.out
    mv t.out $st.raw.assoc.linear
done

metal < p9_metal.sh

# collect the columns together and sort them
# Do not forget to get the allele frequency by running:
# plink --bfile geno --freq

parse_meta_out.py > t.in
head -n 1 t.in > t.out
skh t.in | sort -g -k 14 >> t.out
mv t.out t.in

plot_meta_p_q.R

cp t.in  meta_0.in
cp t.png meta_0.png

# well the first SNP is always AX-83389438
skh t.in | awk '{print $1}' | h -n 1 >> c_snp.ls
# re-run

cp t.in  meta_1.in
cp t.png meta_1.png
# a close match between AX-31642001 and AX-39912161
# because of higher MAF and slightly smaller P
# use AX-31642001
skh t.in | awk '{print $1}' | h -n 1 >> c_snp.ls

# re-run it

cp t.in  meta_2.in
cp t.png meta_2.png

# clearly AX-11541856
skh t.in | awk '{print $1}' | h -n 1 >> c_snp.ls
# re-run it

cp t.in  meta_3.in
cp t.png meta_3.png
# might add AX-31642169
skh t.in | awk '{print $1}' | h -n 1 >> c_snp.ls

# Now the c_snp.ls :
# AX-83389438
# AX-31642001
# AX-11541856
# AX-31642169

# P of the top SNP rise to about 1% afterwards

################################################################################
# 5. Feature selection                                                         #
################################################################################

################################################################################
# 5.1 get SNP dosages

# so far 6 SNPs to be choosen from:
# AX-83389438
# AX-31642001
# AX-11541856
# AX-31642169
# AX-39912161
# AX-11576926

# put them into snp.ls
# or simply the top 20 SNPs in meta 1
skh meta_0.in | h -n 20 | awk '{print $1}' > snp.ls

plink --bfile geno --recode12 --extract snp.ls
get_dosage.py > p9.dosage

################################################################################
# 5.2 select SNP hits

select_SNP.R
# The rfe function of the caret package runs forever!
# This is a simple problem. Don't understand why it takes so loooooong.
# It picks only the top SNP AX-83389438, nothing else.

# according to stepAIC
# final models given the 6 SNPs
# direct ldl_c ~ rc + sex + age + stratum +
#                AX-83389438 + AX-39912161 + AX-31642169 + AX-11541856
#
#   rint_ldl_c ~ AX-83389438 + AX-39912161 + AX-31642001 + AX-31642169 +
#                AX-11541856

# final model given the top 20 SNPs
# direct ldl_c ~ rc + sex + age + stratum +
#                AX-83389438 + AX-39912161 + AX-31642169 + AX-11541856 +
#                AX-50958112 + AX-11447635
#   rint_ldl_c ~ AX-83389438 + AX-39912161 + AX-31642001 + AX-31642169 +
#                AX-11541856 + AX-50958112 + AX-11447635

################################################################################