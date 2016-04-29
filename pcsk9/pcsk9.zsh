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

# stage4.fam ? where is it?
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
