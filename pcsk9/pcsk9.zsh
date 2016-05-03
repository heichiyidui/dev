################################################################################
#                Finding variants in PCSK9 affecting LDL-C                     #
################################################################################

################################################################################
# The PCSK9 gene blocks LDL (Low-density lipoprotein) receptors. Less LDL      #
# receptor on the surface of liver to remove LDL from bloodstream, leads to    #
# higher LDL cholesterol (LDL-C) concentrations.                               #
# PCSK9 is a fat-control drug target. It is well studied with caucasian        #
# populations.                                                                 #
# Now verify our CKB dataset on the detection of SNPs associating with LDL-C   #
# level in this gene.                                                          #
################################################################################


################################################################################
# 1. genotype data PCA

#######################################
# 1.1 the original set

# 32435 subjects in GWAS_sample_ascertainment.txt
# let's start from the stage3 set at
/kuser/shared/data/GWASphase12

tail -n +2 GWAS_sample_ascertainment.txt | awk '{print $2}' > t.ls
grab -f t.ls -c 2 /kuser/shared/data/GWASphase12/stage3.fam | \
    awk '{print $2}' | sort | uniq -d
# subject CK30556316 repeated
# according to
# /kuser/shared/data/GWASphase12/stage3_mandatory_exclusions.txt
# 'NOR15122311 CK30556316 0 0 2 -9' is to be removed

tail -n +2 GWAS_sample_ascertainment.txt | awk '{print $2}' > t.ls
sort t.ls | uniq -d
# no repeat in ck_id

grab -f t.ls -c 2 /kuser/shared/data/GWASphase12/stage3.fam | \
    grep -v "NOR15122311 CK30556316" > t.fam
# 32435 subjects

plink --bfile /kuser/shared/data/GWASphase12/stage3 \
      --keep t.fam --make-bed --out ckb_ph12_s3_qc1
# 659231 SNPs, 32205 people
# CK28710771 got 275 275 het. haploid genotypes present

#######################################
# 1.2 LD based pruning

plink  --bfile ckb_ph12_s3_qc1 \
    --autosome \
    --geno 0.01 \
    --hwe 1e-4 midp \
    --maf 0.05 \
    --indep-pairwise 1500 150 0.2 \
    --make-bed --out ckb_ph12_s3_qc1_ldfree

# 636670 SNPs loaded(autosome), 32435 subjects
#  73412 variants removed due to missing genotype data (--geno).
#  17297 variants removed due to Hardy-Weinberg exact test (--hwe).
# 208608 variants removed due to minor allele threshold(s)

# 337353 SNPs in, 217191 pruned, 120162 left

plink --bfile ckb_ph12_s3_qc1_ldfree \
    --extract ckb_ph12_s3_qc1_ldfree.prune.in \
    --make-bed --out  ckb_ph12_s3_qc1_ldfree
# 120162 SNPs, 32435 people

#######################################
# 1.3 First PCA

tail -n +2 GWAS_sample_ascertainment.txt | awk '{print $3}' | sort | uniq -d
# no repeat in study ids
tail -n +2 GWAS_sample_ascertainment.txt | awk '{print substr($3,1,2)}' | \
    sort | uniq -c
# 1388 12 Qingdao
# 3310 16 Harbin
# 1174 26 Haikou
# 1705 36 Suzhou
# 2459 46 Liuzhou
# 4106 52 Sichuan
# 4705 58 Gansu
# 4343 68 Henan
# 3363 78 Zhejiang
# 5882 88 Hunan
# all from the 10 regions





random_shuffle_lines.py ckb_ph12_s3_qc1_ldfree.fam | h -n 100 > t.fam

plink --bfile ckb_ph12_s3_qc1_ldfree --keep t.fam --make-bed --out t100



/kuser/shared/bin/EIG/bin/smartpca.perl \
    -i t100.bed \
    -a t100.bim \
    -b t100.fam \
    -o t100_o   \
    -p t100_o.plot \
    -e t100_o.eval \
    -l t100_o.log

601.57s user 65.66s system 180% cpu 6:10.50 total
# 6 min for 2000 individuals
# about 5 hours on 32435 individuals

# GSL is missing on the NC2! Emailed the administrator.
# Mike Weale's script of running EIGENSOFT
# EIGENSOFTplus_v12.r
wget https://raw.githubusercontent.com/KHP-Informatics/chip_gt/master/\
protocols_and_QC/exome_chip_QC/EIGENSOFTplus_v12.r

R --vanilla --slave --args stem=t100 \
    ESOFTdir=/kuser/shared/bin/EIG \
    < EIGENSOFTplus_v12.r

R --vanilla --slave --args stem=ckb_ph12_s3_qc1_ldfree \
    ESOFTdir=/home/kuang/bin/EIG \
    < EIGENSOFTplus_v12.r

# get a random set of 100 subjects
plink --bfile t100 --pca
