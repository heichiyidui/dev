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

################################################################################
# 1. genotype data                                                             #
################################################################################

################################################################################
# 1.1 the original set

# Start from the stage3 set at
/kuser/shared/data/GWASphase12

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

# Do not perform it yet.

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

# Do not remove subjects yet

plink --bfile ckb_ph12_s3_qc01 \
      --het

# No, plink says use a LD-free set.
# I made a comparison. The two sets of F valuse are pretty close.
# But let's just do the LD-pruning

plink --bfile ckb_ph12_s3_qc01 \
      --geno 0.01 \
      --hwe 1e-4 midp \
      --maf 0.05 \
      --indep-pairwise 1500 150 0.2

# 120201 SNPs left

plink --bfile ckb_ph12_s3_qc01 \
      --extract plink.prune.in \
      --make-bed --out t

plink --bfile t  \
      --het

# and missing
plink --bfile ckb_ph12_s3_qc01 \
      --missing

paste plink.het plink.imiss | awk '{print $1,$2,$6,$12}' > t.in

#--------------------------------------
# in R

library(ggplot2)
data = read.table('t.in',header=T)

p1 <- ggplot(data) + geom_point(aes(x=F, y=F_MISS)) + theme_bw()

p2 <- ggplot(data=subset(data,F<0.10), aes(F)) +
      geom_histogram(bins=30) + theme_bw()

threshold = mean(data$F) - 3 * sd(data$F)
data$lowF <- data$F < threshold
# -0.0287 ~ 0.03048 for +-3 SD
# 395 individuals flaged
# or 20 individuals on the lower end

p3 <- ggplot(data,aes(x=lowF, y=F_MISS)) + geom_boxplot( ) + theme_bw()

write.table(subset(data,lowF)$IID,file='low_hom.ls',
            row.names=F,col.names=F,quote=F)
#--------------------------------------
# end of the R script

grab -f low_hom.ls -c 2 ckb_ph12_s3_qc01.fam | awk '{print $1}' | sort | uniq -c
# The plate ids are rather scattered. I was expecting some plate contamination.

# low homozygotes might mean sample contamination.
# The 20 subjects with < 3 x SD F values are of very high missingness.
# 0.022 vs 0.0029 for all subjects
# And, according to the plink IBD estimation,
# these subjects are related to almost EVERYONE else.
# templed to remove them

#######################################
# don't remove any subject here yet

# grab -f low_hom.ls -c 2 ckb_ph12_s3_qc01.fam > t.fam

# plink --bfile ckb_ph12_s3_qc01 \
#       --remove t.fam \
#       --make-bed --out  ckb_ph12_s3_qc02

# # 553726 variants and 32185 people left

################################################################################
# 1.4 IBD and PCA

#######################################
# 1.4.1 LD pruning again

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

# 20 minutes on NC2
# plink.genome is huge without the 0.05 filter

tail -n +2 plink.genome | awk '{print $2 "\n" $4}'  | sort | uniq > t.ls
# 32040 subjects to be removed? not acceptable.

tail -n +2 plink.genome | awk '{print $2,$4}' > t.in

select_fam.py > to_remove.ls
# Remove (flag) the subjects according to their degree of connection.
# The most related subjects are to be removed first.
# The top ones are still with low homozygotes.

# 6990 subjects, 22% of 32205

# There are three pairs of people with PI_HAT about 0.5
# CK22775935 CK28153131
# CK28568808 CK28569055
# CK30580964 CK30589772
# The others are a huge family.

# We need to delete all edges. I tried different ad hoc approaches to remove as
# few as possible nodes. I Was doing the PCA with 7333 subjects flagged as
# related. It shouldn't change the results much.

# change the 6th column of pca.fam to 'rel' and 'no_rel'.
grab -f to_remove.ls -c 2 pca.fam | awk '{print $1,$2,$3,$4,$5,"rel"}' > t.out
grab -f to_remove.ls -c 2 pca.fam -v | \
    awk '{print $1,$2,$3,$4,$5,"no_rel"}' >> t.out

awk '{print $2}' pca.fam > t.ls
sort_table -f t.ls -c 2 t.out > t2.out

mv t2.out pca.fam

printf "no_rel" > pop.ls

#######################################
# 1.4.3 PCA via EIGENSOFT

# By default, smartpca should be using multithreading now.
# Runing PCA with the un-related subjects only. Then project the factors to the
# related subjects as well.
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

awk '{print $1,$2,$6}' pca.fam > pca_cluster.in

nohup plink --bfile pca \
            --pca 10 \
            --within pca_cluster.in \
            --pca-cluster-names no_rel &

# plink uses about 80 threads
# 7G of memory
# about 5 hours

################################################################################
# 1.5 final genotype QC

#######################################
# 1.5.1 check the PCA plots

sed 's/\:/\t/' no_rel.pca.evec -i

printf "CK24820387\nCK25228869\nCK28902540\nCK28730586" > t.ls

tail -n +2 no_rel.pca.evec |\
    awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | \
    grep -f t.ls -v > t.in

# for ck_id.ls and study_id.ls, see section 2.
sort_table -f ck_id.ls t.in > t.out

echo "study_id rc pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 is_fam" > t.in

paste study_id.ls t.out | \
    awk '{print $1,substr($1,1,2),$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' >> t.in

plot_pca.R
cp t.in pca.in

# It seems no outliers are left to be removed.
# PCs are not corresponding to if the subjects are related or not.
# PCs, up to PC7 and PC8, are obviously related to the region codes.

# If we didn't remove the Suzhou (and other) families, they will cousing trouble
# with PCA. PC4-PC10 are all draged by them.


# plink is about 4 times faster than eigensoft.
# They produce mostly the same results.
#

#######################################
# 1.5.2 Remove badly called SNPs found in the manual check.

# Don't do it yet.

# tail -n +2 ../ckb_batch_check/manual_chk_res.table | \
#     awk '{if ($2==0) print $1}' > to_remove_snp.ls

# # 13938 SNPs to be removed.
# # and the four subjects with missing ascertainments

# printf "CK24820387\nCK25228869\nCK28902540\nCK28730586" > t.ls

# grab -f t.ls -c 2 ckb_ph12_s3_qc02.fam > t.fam

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

printf "CK24820387\nCK25228869\nCK28902540\nCK28730586" > t.ls
grep -f t.ls -v GWAS_SNPdata_samples.csv > t.out
mv t.out GWAS_SNPdata_samples.csv

tail -n +2  GWAS_SNPdata_samples.csv | awk -F"," '{print $1}' > ck_id.ls
tail -n +2  GWAS_SNPdata_samples.csv | awk -F"," '{print $2}' > study_id.ls

# 32201 uniq study and ck id pairs.
# ck_id.ls and study_id.ls were used in 1.5.1

################################################################################
# 2.2 stratification

# use the sheet 2 of PCSK9_sample_summary.xlsx from
# K:\kadoorie\Groups\Genetics\PROJECTS\PCSK9
# for subject stratification.
# Added 1288 '0's in the 'pass GWAS QC' column according to study_id.ls
# instead of 1275 in the original one.

# simplified table with header:
# studyid ascert still_OK dir_ldl_base dir_ldl_rs1 dir_ldl_rs2 indir_ldl_rs2

# use LDL-c_biochem_data.xlsx for direct LDL-C
# use ldl_levels_resurvey2_latest.xls for indirect LDL-C
# scale indirect ldl-c by 0.01
# 18091 direct and 7110 indirect measures

# indirect ldl-c is smaller
# mean 2.120 vs 2.376
# t-test says the difference is very significant.

PCSK9_sample_summary.csv

# put the first 7 columns of PCSK9_sample_summary.csv into t.in
# and prepare the other input files:
# age_base.csv
# age_resu1.csv
# age_resu2.csv
# direct_ldl_c.txt
# indirect_ldl_c.txt

get_strat.py > t.out


# str1: ICH     direct       4762
# str2: IS      direct       5210
# str3: SAH     direct        167
# str4: MI/IHD  direct       1265
# str5: control direct       6696
# str6: all indirect left    4174

# others: 11205 NA


################################################################################
# 2.3  Covariates and phenotypes

#######################################
# 2.3.1 The phenotype file is easy to make.
pheno.csv
# FID IID LDL
# 21553 subject LDL measures

#######################################
# 2.3.2 The covariates file

# use plink no-rel 10 PCs

# using the data request form to get the age_at_study info.
# ages of 3 subjects are missing:
# CK28728060      580282304
# CK28728462      580281490
# CK30579300-1    580235861

# 1377 12 Qingdao   rc1
# 3285 16 Harbin    rc2
# 1164 26 Haikou    rc3
# 1695 36 Suzhou    rc4
# 2424 46 Liuzhou   rc5
# 4090 52 Sichuan   rc6
# 4675 58 Gansu     rc7
# 4324 68 Henan     rc8
# 3336 78 Zhejiang  rc9
# 5831 88 Hunan     Hunan will be absent from RC factors

printf "CK24820387\nCK25228869\nCK28902540\nCK28730586" > t.ls

tail -n +2 no_rel.pca.evec |\
    awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | \
    grep -f t.ls -v > t.in

paste study_id.ls t.in | sed 's/\t/\ /' > t.out

# ...

cov.csv

# 'NA' to -9 for missing covariates (age)
# Given the stratification only 4 subjects have missing ages.

# don't forget keep-pheno-on-missing-cov

################################################################################
# 3. plink linear regression                                                   #
################################################################################

# get the sets
awk '{if ($8==1) print $1}' PCSK9_sample_summary.csv > t1.ls ;
grab -f t1.ls -c 2 pca.fam > t1.fam ;
plink --bfile ckb_ph12_s3 --keep t1.fam --make-bed --out st1 &
awk '{if ($8==2) print $1}' PCSK9_sample_summary.csv > t2.ls ;
grab -f t2.ls -c 2 pca.fam > t2.fam ;
plink --bfile ckb_ph12_s3 --keep t2.fam --make-bed --out st2 &
awk '{if ($8==3) print $1}' PCSK9_sample_summary.csv > t3.ls ;
grab -f t3.ls -c 2 pca.fam > t3.fam ;
plink --bfile ckb_ph12_s3 --keep t3.fam --make-bed --out st3 &
awk '{if ($8==4) print $1}' PCSK9_sample_summary.csv > t4.ls ;
grab -f t4.ls -c 2 pca.fam > t4.fam ;
plink --bfile ckb_ph12_s3 --keep t4.fam --make-bed --out st4 &
awk '{if ($8==5) print $1}' PCSK9_sample_summary.csv > t5.ls ;
grab -f t5.ls -c 2 pca.fam > t5.fam ;
plink --bfile ckb_ph12_s3 --keep t5.fam --make-bed --out st5 &
awk '{if ($8==6) print $1}' PCSK9_sample_summary.csv > t6.ls ;
grab -f t6.ls -c 2 pca.fam > t6.fam ;
plink --bfile ckb_ph12_s3 --keep t6.fam --make-bed --out st6 &

# linear regression with all covariants
for st in st1 st2 st4 st5 st6 ; do
    nohup plink --bfile $st \
      --pheno pheno.csv \
      --pheno-name LDL \
      --covar cov.csv keep-pheno-on-missing-cov \
      --linear hide-covar --ci 0.95 \
      --out $st &
done

# The SAH cohort is way too small.
# With all covariants the output will be all NA.
nohup plink --bfile st3 \
  --pheno pheno.csv \
  --pheno-name LDL \
  --covar cov.csv keep-pheno-on-missing-cov \
  --covar-name sex,age,pc1,pc2 \
  --linear hide-covar --ci 0.95 \
  --out st3 &

#######################################
# format the output

for st in st1 st2 st3 st4 st5 st6 ; do
    head -n 1 $st.assoc.linear > t.out
    grep ADD $st.assoc.linear | grep -v NA >> t.out
    mv t.out $st.assoc.linear
done

# put A2 into the tables
for st in st1 st2 st3 st4 st5 st6 ; do
    add_a2.py $st.assoc.linear > t.out
    mv t.out $st.assoc.linear
done

#######################################
# QQ and Manhattan plots and lambda

# let's have a look
for st in st1 st2 st3 st4 st5 st6 ; do
    plot_qq_man.R $st.assoc.linear &
done

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

for st in st1 st2 st3 st4 st5 st6 ; do
    lambda.R $st.assoc.linear
done

# st1.assoc.linear lambda:  1.007485
# st2.assoc.linear lambda:  1.008424
# st3.assoc.linear lambda:  0.9920914
# st4.assoc.linear lambda:  0.9907005
# st5.assoc.linear lambda:  1.025426
# st6.assoc.linear lambda:  1.008893

################################################################################
# 4. METAL analysis                                                            #
################################################################################

plink --meta-analysis  st6.assoc.linear st3.assoc.linear + qt



# prepare for plotting

awk '{print $1}' pcsk9_direct1.tbl | t -n +2 > t.ls
grab -f t.ls -c 2 ckb_ph12_s3.bim | awk '{print $1,$2,$4} ' > t.in
sort_table -f t.ls -c 2 t.in > t.out

mv t.out t.in
t -n +2 pcsk9_direct1.tbl > t2.in
printf "CHR SNP BP A1 A2 BETA SE P DIR" > pcsk9_direct_metal.out
paste t.in t2.in | \
     awk '{print $1,$2,$3,$5,$6,$7,$8,$9,$10}'>> pcsk9_direct_metal.out




awk '{print $1}' pcsk9_all1.tbl | t -n +2 > t.ls
grab -f t.ls -c 2 ckb_ph12_s3.bim | awk '{print $1,$2,$4} ' > t.in
sort_table -f t.ls -c 2 t.in > t.out

mv t.out t.in
t -n +2 pcsk9_all1.tbl > t2.in
printf "CHR SNP BP A1 A2 BETA SE P DIR" > pcsk9_all_metal.out
paste t.in t2.in | \
     awk '{print $1,$2,$3,$5,$6,$7,$8,$9,$10}'>> pcsk9_all_metal.out







################################################################################





#######################################
# 1.3 LD-based pruning

plink  --bfile ckb_ph12_s3_qc02 \
    --geno 0.01 \
    --hwe 1e-4 midp \
    --maf 0.05 \
    --indep-pairwise 1500 150 0.2

# 522846 SNPs loaded, 32109 subjects
#  15412 SNPs removed due to missing genotype data (--geno).
#   7298 SNPs removed due to Hardy-Weinberg exact test (--hwe).
# 164073 SNPs removed due to minor allele threshold(s)

# 336063 SNPs in, 216173 pruned, 119890 left

plink --bfile ckb_ph12_s3_qc02 \
    --extract plink.prune.in \
    --make-bed --out  ckb_ph12_pca
# 119890 SNPs, 32435 people

#######################################
# 1.4 split the plink files into 10 RC sets

for rc in 12 16 26 36 46 52 58 68 78 88 ; do
    tail -n +2 GWAS_sample_ascertainment.txt | \
        awk -v rc="$rc" '{if (substr($2,1,2)==rc) print $1}' > rc_id.ls
    grab -f rc_id.ls -c 2 ckb_ph12_pca.fam > t.fam
    plink --bfile ckb_ph12_pca \
          --keep t.fam \
          --make-bed --out ckb_ph12_pca_rc$rc
done

# Pretend all are control for the moment
for rc in 12 16 26 36 46 52 58 68 78 88 ; do
    awk '{$6=1; print $0}' ckb_ph12_pca_rc$rc.fam > t.fam
    mv t.fam  ckb_ph12_pca_rc$rc.fam
done

#######################################
# 1.5 families clusters


for rc in 12 16 26 36 46 52 58 68 78 88 ; do
    nohup plink --bfile ckb_ph12_pca_rc$rc \
          --genome --out rc$rc  &
done

# according to genome ibs distances,
# there are relatives in the cohorts, but no duplications now.

for rc in 12 16 26 36 46 52 58 68 78 88 ; do
    echo $rc
    tail -n +2 rc$rc.genome | awk '{if ($10>0.40) print $2 "\n" $4}' | \
        sort | uniq | wc
done

# some families to be removed during PCA
# using the threshold phat >= 0.125 for adding relatives

for ifile in *.genome ; do
    echo $ifile
    count_fam.py $ifile | sort -g | tail -n 3
done

# 12 Qingdao  1388      4 in the lastest family
# 16 Harbin   3310      3 in the lastest family
# 26 Haikou   1173     10 in the lastest family
# 36 Suzhou   1705    129 in the lastest family
# 46 Liuzhou  2445      4 in the lastest family
# 52 Sichuan  4106     12 in the lastest family
# 58 Gansu    4703     11 in the lastest family
# 68 Henan    4341      9 in the lastest family
# 78 Zhejiang 3360     18 in the lastest family
# 88 Hunan    5881     12 in the lastest family

# put the 129 ids of the Suzhou family into 'suzhou_fam.ls'
# edit the last column of ckb_ph12_pca_rc36.fam
# so that we have 'fam_1' or 'no_rel' in the last column

echo "no_rel" > t.pop
# population list file for smartpca on rc36

# for rc in 12 16 26 36 46 52 58 68 78 88 ; do
#     plink --bfile ckb_ph12_pca_rc$rc \
#           --read-genome rc$rc.genome \
#           --cluster \
#           --out rc$rc
# done
# all single clusters, mostly useless
# Should be some way to tweak it.

#######################################
# 1.6 first round region PCA

# no rc36
for rc in 12 16 26 46 52 58 68 78 88 ; do
    ~/bin/EIG/bin/smartpca.perl \
        -i ckb_ph12_pca_rc$rc.bed \
        -a ckb_ph12_pca_rc$rc.bim \
        -b ckb_ph12_pca_rc$rc.fam \
        -o rc$rc.pca  \
        -p rc$rc.plot \
        -e rc$rc.eval \
        -l rc$rc.log
done

for rc in 36 ; do
    ~/bin/EIG/bin/smartpca.perl \
        -i ckb_ph12_pca_rc$rc.bed \
        -a ckb_ph12_pca_rc$rc.bim \
        -b ckb_ph12_pca_rc$rc.fam \
        -w t.pop \
        -o rc$rc.pca  \
        -p rc$rc.plot \
        -e rc$rc.eval \
        -l rc$rc.log
done

# 6 min for 2000 individuals
# about 5 hours on 32435 individuals
# But out of memory on a PC with 12G memory.

# The GSL library is missing on the NC2! Emailed the administrator.

# So far only rc 46 PCA shows some internal structure there.
# Use the first 10 global PCs, from stage4_pca.xlsx

ckb_global_pca.txt

################################################################################
# 2 LDL-C measurements

# from
# file://K:\kadoorie\Groups\CKB-Statistics\Biochemistry Data\20160419_final_data
# 20160426_101958_biochemistry_result_grid.xls

# about 20K subjects
# removed the 'diluted'
# removed the 'Result may be affected by lipemia' etc

pheno.ods
# we need their DOB.
# we need to add sttroke and statin use.

# checking different covariants
plot_pheno.R

# put the ldl phenotype into ldl.pheno

# rc cov as 9 binary factors
# rc_map = {'12':0, '16':1, '26':2, '36':3, '46':4,
#           '52':5, '58':6, '68':7, '78':8}

# ascertainment cov as 4 binary factors
# rc_map = {'control':0, 'ICH':1, 'IS':2, 'MI_IHD':3}

# put the covariants into ldl.cov

################################################################################
# 3 genotype file

# almost without any QC

plink --bfile ckb_ph12_s3_qc01 \
      --chr 1 --from-mb 55.44 --to-mb 56.01 \
      --keep ldl.pheno \
      --make-bed --out geno
# 159 SNPs, 17458 subjects


################################################################################
# 4. assoc

plink --bfile geno \
      --pheno ldl.pheno \
      --linear --ci 0.95

# 1   AX-83389438   55509585    T        ADD
# 17455   -0.09347 0.007536  -0.1082 -0.07869        -12.4    3.588e-35

# with all covariants
plink --bfile geno \
      --pheno ldl.pheno \
      --covar ldl.cov \
      --linear --ci 0.95

grep ADD plink.assoc.linear | sort -g -k 12 | head -n 22
head  -n 1 plink.assoc.linear > t.in
grep ADD plink.assoc.linear | sort -g -k 12 | grep -v NA  >> t.in
# then plot it

#   1   AX-83389438   55509585    T        ADD
# 17455      -0.38  0.03056  -0.4399  -0.3201       -12.43    2.419e-35
# MAF = 0.013


# check ld
plink --bfile geno \
       --r2 --ld-snp AX-83389438 \
       --ld-window-r2 0 \
       --ld-window 99999 \
       --ld-window-kb 77000
# nothing is in LD with AX-83389438

sort -g -k 7 plink.ld
#   1     55509585   AX-83389438      1     55601339   AX-11187195    0.0140177
# highest of all but itself is 0.0140177

plink --bfile geno \
       --r2 --ld-snp AX-11576926 \
       --ld-window-r2 0 \
       --ld-window 99999 \
       --ld-window-kb 77000
# not too bad

plink --bfile geno \
       --r2 --ld-snp AX-39912161 \
       --ld-window-r2 0 \
       --ld-window 99999 \
       --ld-window-kb 77000
#############
# chk the clustering
