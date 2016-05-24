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

# Obviously, the gender check has been done.

################################################################################
# 1.2 very basic missing-call and MAF filter

plink --bfile ckb_ph12_s3 \
      --geno 0.05 \
      --maf  0.0001 \
      --make-bed --out ckb_ph12_s3_qc01

# 51983 variants removed due to missing genotype data (--geno).
# 30961 variants removed due to minor allele threshold(s)
# 553726 variants and 32205 people left

################################################################################
# 1.3 heterogeneous related to missingness

plink --bfile ckb_ph12_s3_qc01 \
      --het

# No, plink says use a LD-free set.
# I made a comparison, the two sets of F valuse are pretty close.
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
R
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

grab -f low_hom.ls -c 2 ckb_ph12_s3_qc01.fam | awk '{print $1}' | sort | uniq -c
# The plate ids are rather scattered. I was expecting some plate  contamination.

#--------------------------------------
#

# low homozygotes might mean sample contamination.
# The 20 subjects with < 3 x SD F values are of very high missingness.
# 0.022146 vs 0.002926294 for all
# And, according to the plink IBD estimation,
# these subjects are related to EVERYONE else.

grab -f low_hom.ls -c 2 ckb_ph12_s3_qc01.fam > t.fam

plink --bfile ckb_ph12_s3_qc01 \
      --remove t.fam \
      --make-bed --out  ckb_ph12_s3_qc02

# 553726 variants and 32185 people left

################################################################################
# 1.4 IBD and PCA

# LD pruning again

plink --bfile ckb_ph12_s3_qc02 \
      --geno 0.01 \
      --hwe 1e-4 midp \
      --maf 0.05 \
      --indep-pairwise 1500 150 0.2
# 120177 SNPs in, 217508 SNPs out

plink --bfile ckb_ph12_s3_qc02 \
      --extract plink.prune.in \
      --make-bed --out pca
# 120177 variants and 32185 people

plink --bfile pca \
      --genome

# 20 minutes on NC2

awk '{if ($10> 0.05) print $0}' plink.genome > t_genome.in
rm plink.genome
# plink.genome is huge
# now t_genome.in is much much smaller

awk '{if ($10>0.05) print $2 "\n" $4}' t_genome.in | sort | uniq > t.ls
# 31771 subjects to be removed? not acceptable.

awk '{if ($10>0.05) print $2 "\n" $4}' t_genome.in | sort | uniq -c | \
    sort -k 1 -g -r > t.dat
# sort the subject ids according to the frequencies they are related to others.
# The most related subjects are still with low homozygotes.

select_fam.py > to_remove.ls
# select ck_ids to be removed. The most related subjects to be removed first.
# 7191 subjects, 22.3% of 32185

grab -f to_remove.ls -c 2 pca.fam | awk '{print $1,$2,$3,$4,$5,"rel"}' > t.out
grab -f to_remove.ls -c 2 pca.fam -v | \
    awk '{print $1,$2,$3,$4,$5,"no_rel"}' >> t.out

awk '{print $2}' pca.fam > t.ls
sort_table -f t.ls -c 2 t.out > t2.out

mv t2.out pca.fam
# change the 6th column of pca.fam to 'rel' and 'no_rel'.

echo "no_rel" > pop.ls

# By default, smartpca should be using multithreading now.
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

nohup /kuser/shared/bin/EIG/bin/smartpca.perl \
        -i pca.bed \
        -a pca.bim \
        -b pca.fam \
        -o all.pca  \
        -p all.plot \
        -e all.eval \
        -l all.log  \
        -m 0   &

################################################################################
# 1.5

# check the PCA plots and don't forget to remove the bad SNPs found via
# the manual checking of the plate and batch effects.


################################################################################
# 2. phenotype data                                                            #
################################################################################

# 32435 subjects in GWAS_sample_ascertainment.txt from
# K:\kadoorie\Groups\Genetics\Data Archive\Project Sample Lists\Lists\
# GWAS_SNPdata_samples.xlsx

# I replaced ' ' with '_' in the 'ascert.' and 'notes' columns.


awk '{print $2}' ckb_ph12_s3_qc01.fam > ck_id.ls
# 32109 uniq ids in ck_id.ls


# direct LDL-C
# min:       0.3300
# max:       9.9400
# mean:      2.3782
# std:       0.7136
# quartiles: 0.3300 1.8800 2.3200 2.8000 9.9400

# indirect LDL-C
# min:       0.0000
# max:       1078.0000
# mean:      218.6738
# std:       91.0916
# quartiles: 0.0000 153.0000 212.0000 274.0000 1078.0000

# just use indirect LDL-C / 100

awk '{if ($2=="ICH" && $3==0 && $5==1) print $1}' t.in > st1.ls
awk '{if ($2=="IS"  && $3==0 && $5==1) print $1}' t.in > st2.ls
awk '{if ($2=="SAH" && $3==0 && $5==1) print $1}' t.in > st3.ls
awk '{if ($2=="MI/IHD" && $3==0 && $5==1) print $1}' t.in > st4.ls
awk '{if ($2=="control" && $3==0 && $5==1) print $1}' t.in > st5.ls
awk '{if ($2=="resurvey_2" && $3==0 && $5!=1 && $8 != "NA" ) print $1}' \
  t.in > st6.ls

wc st?.ls
#   4694   4694  51998 st1.ls
#   5113   5113  56683 st2.ls
#    165    165   1823 st3.ls
#   1253   1253  13965 st4.ls
#   5756   5756  63715 st5.ls
#   3201   3201  35833 st6.ls
#  20182  20182 224017 total


# 10 RCs
# 1388 12 Qingdao
# 3310 16 Harbin
# 1173 26 Haikou
# 1705 36 Suzhou
# 2445 46 Liuzhou
# 4106 52 Sichuan
# 4703 58 Gansu
# 4341 68 Henan
# 3360 78 Zhejiang
# 5881 88 Hunan

# 10 PCs

# all in cov.csv


################################################################################
# 3. linear regression


# with all covariants
plink --bfile geno \
      --pheno ldl.pheno \
      --covar ldl.cov \
      --linear --ci 0.95

grab -f st1.ls ckb_ph12_s3_qc01.fam -c 2 > t1.fam ; \
plink --bfile ckb_ph12_s3_qc01 --keep t1.fam --make-bed --out st1 &
grab -f st2.ls ckb_ph12_s3_qc01.fam -c 2 > t2.fam ; \
plink --bfile ckb_ph12_s3_qc01 --keep t2.fam --make-bed --out st2 &
grab -f st3.ls ckb_ph12_s3_qc01.fam -c 2 > t3.fam ; \
plink --bfile ckb_ph12_s3_qc01 --keep t3.fam --make-bed --out st3 &
grab -f st4.ls ckb_ph12_s3_qc01.fam -c 2 > t4.fam ; \
plink --bfile ckb_ph12_s3_qc01 --keep t4.fam --make-bed --out st4 &
grab -f st5.ls ckb_ph12_s3_qc01.fam -c 2 > t5.fam ; \
plink --bfile ckb_ph12_s3_qc01 --keep t5.fam --make-bed --out st5 &
grab -f st6.ls ckb_ph12_s3_qc01.fam -c 2 > t6.fam ; \
plink --bfile ckb_ph12_s3_qc01 --keep t6.fam --make-bed --out st6 &



# with all covariants
for st in st1 st2 st4 st5 st6 ; do
    nohup plink --bfile $st \
      --pheno ldl.pheno \
      --covar ldl.cov \
      --linear --ci 0.95 \
      --out $st &
done

# SAH is way too small, witch all covariants the output will only be NA
plink --bfile st3 --pheno ldl.pheno --covar ldl.cov \
    --covar-name age, gender, PC1, PC2 \
    --linear --ci 0.95 --out st3 &

for st in st1 st2 st3 st4 st5 st6 ; do
    head -n 1 $st.assoc.linear > $st.out ;
    grep ADD $st.assoc.linear >> $st.out &
done

# lets have a look
for st in st1 st2 st3 st4 st5 st6 ; do
    plot_qq_man.R $st  &
done


################################################################################
# 4. METAL analysis


# put a2 into it

for st in st1 st2 st3 st4 st5 st6 ; do
    t1.py $st > $st.metal_in  &
done


# prepare for plotting

awk '{print $1}' pcsk9_direct1.tbl | t -n +2 > t.ls
grab -f t.ls -c 2 ckb_ph12_s3.bim | awk '{print $1,$2,$4} ' > t.in
sort_table -f t.ls -c 2 t.in > t.out

mv t.out t.in
t -n +2 pcsk9_direct1.tbl > t2.in
echo "CHR SNP BP A1 A2 BETA SE P DIR" > pcsk9_direct_metal.out
paste t.in t2.in | \
     awk '{print $1,$2,$3,$5,$6,$7,$8,$9,$10}'>> pcsk9_direct_metal.out




awk '{print $1}' pcsk9_all1.tbl | t -n +2 > t.ls
grab -f t.ls -c 2 ckb_ph12_s3.bim | awk '{print $1,$2,$4} ' > t.in
sort_table -f t.ls -c 2 t.in > t.out

mv t.out t.in
t -n +2 pcsk9_all1.tbl > t2.in
echo "CHR SNP BP A1 A2 BETA SE P DIR" > pcsk9_all_metal.out
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
