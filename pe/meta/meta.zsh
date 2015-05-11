################################################################################
# After the Wednesday meeting, 
# TODO:
#
#1) Check the accuracy of the pedigrees using the genotyping data and do further
# QC accordingly. Chris to advise on suitable programme/s please.
#
#2) A meta-analysis within each centre in the PE sample to check the PCA
# ancestry covariates we've used. For this we'll use the unrelated case-control 
# dataset.
#
#3) Look at the ancestry PCA graph separating cases from controls.
#
#4) Complete the analysis including the additional 1072 samples from the 
# families using Unphased.
#
#5) Test published SNPs for schizophrenia in our sample. To do this with the 
# broad phenotype "psychosis" and also to look at the subset with the stringent 
# phenotype "schizophrenia". We'll try to do the same with the subset with 
# "non-schizophrenia psychoses" but the sample may be too small. Compare the 
# evidence for association looking at effect sizes (OR, not just the p values).
#
#6) We'll do a preliminary analysis of the endophenotypes.
#
#7) Think about publication strategy: A) Publish case-control analysis first and
# endophenotypes as a second paper. B) Produce one paper including both 
# case-control and endophenotypes.
################################################################################

# 1. meta-analysis of cohorts

# There are 8 cohorts in our set
# PE_EDIN     ed
# PE_HEID     he
# PE_HOLL     ho
# PE_LOND     lo
# PE_MUNI     mu
# PE_PAMP     pa
# Valdecilla  va
# WAFSS       wa 

# mu has only controls
# pa has only patients

# combine mu with he
# combine pa with va

# pe17_ed
# pe17_muhe
# pe17_ho
# pe17_lo
# pe17_pava
# pe17_wa

plink --noweb --bfile pe17_ed   --out pe17_ed   --logistic --ci 0.95 --adjust & 
plink --noweb --bfile pe17_muhe --out pe17_muhe --logistic --ci 0.95 --adjust & 
plink --noweb --bfile pe17_ho   --out pe17_ho   --logistic --ci 0.95 --adjust & 
plink --noweb --bfile pe17_lo   --out pe17_lo   --logistic --ci 0.95 --adjust & 
plink --noweb --bfile pe17_pava --out pe17_pava --logistic --ci 0.95 --adjust & 
plink --noweb --bfile pe17_wa   --out pe17_wa   --logistic --ci 0.95 --adjust & 

plink  --noweb --meta-analysis \
pe17_ed.assoc.logistic  pe17_lo.assoc.logistic    pe17_pava.assoc.logistic \
pe17_ho.assoc.logistic  pe17_muhe.assoc.logistic  pe17_wa.assoc.logistic \
--out meta

t1.R
----------------------------------------
source("http://www.StephenTurner.us/qqman.r")

results <- read.table("meta.meta",T);
png(filename='meta_man.png',width=1600,height=600);manhattan(results);
dev.off();
png('meta_qq.png',width=1200,height=1200);qq(results$P);dev.off();
----------------------------------------
plink --noweb --bfile ../assoc/pe17 --logistic --out pe17
plink --noweb --bfile ../assoc/pe17 --logistic --covar ../assoc/t.cov \
 --covar-number 1-3 --out pe17_cov3 --adjust

# cov0  : 695464 SNPs, 448 with less than 1^-4 p values 
# cov3  : 695464 SNPs, 114 with less than 1^-4 p values
# meta6 : 683557 SNPs, 84  with less than 1^-4 p values
# upcv3 : 695464 SNPs, 61  with less than 1^-4 p values

# cov3 (74 ((40)) 44) meta6
# 
# upc3 (53 ((8 )) 76) meta6
# 
# upc3 (40 ((21)) 93) cov3

# cov3 (44 ((70)) 378) cov0
# 
# upc3 (40 ((21)) 427) cov0
# 
# meta6(44 ((42)) 406) cov0
