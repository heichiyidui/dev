################################################################################
# Linear regression of the batch 3 data set                                    #
################################################################################

# Load libraries
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(limma)
library(lumiHumanIDMapping)

source("code/sjnewhouse_misc_R.R")

#
load(
'ba3/BA3_lumi_processing_t_out/eset_bg_log2_rsn/BA3.eset_bg_log2_rsn_adj.RData')

eset_bg_log2_rsn_adj

# using Subject_Demographics_with_chip_data_for_processing_April2015_FINAL
# column I and S 
# AD, MCI and CTL only
# 449 samples 
# 57 samples to be removed

samples_to_remove=c(
"9370786090_B","9370786090_D","9370786090_E","9370786091_D","9370786091_H",
"9370786098_B","9370786098_G","9371242013_A","9371242063_C","9402444013_A",
"9402444013_D","9402444013_F","9402444013_H","9402444013_K","9402444016_E",
"9402444017_A","9402444017_C","9402444022_D","9402444022_J","9402444024_C",
"9402444024_L","9402444029_G","9402444029_J","9402444029_K","9402444030_D",
"9402444034_B","9402444040_D","9464921149_B","9477874003_F","9477874019_G",
"9477874020_J","9477874036_B","9477874036_D","9477874050_B","9477874050_J",
"9477874066_D","9477874066_F","9477874066_I","9477874069_I","9477874075_F",
"9477874075_L","9477874079_B","9477874081_C","9477874081_D","9477874081_K",
"9534190031_L","9534190041_G","9534190041_H","9534190041_L","9534190085_C",
"9534190085_H","9534190085_L","9534190110_D","9534190110_F","9534190110_H",
"9534190112_L","9703789038_A")

eset=removeSamples_eset_lumi(eset_bg_log2_rsn_adj,samples_to_remove)
# 449 samples, 3976 probes

pData(eset)$GROUPS = c(
 "AD" ,"AD" ,"MCI","MCI","AD" ,"CTL","MCI","MCI","AD" ,"AD" ,"AD" ,"AD" ,
 "MCI","AD" ,"MCI","MCI","MCI","AD" ,"AD" ,"CTL","AD" ,"CTL","AD" ,"AD" ,
 "MCI","MCI","AD" ,"AD" ,"AD" ,"AD" ,"CTL","AD" ,"AD" ,"AD" ,"AD" ,"CTL",
 "CTL","AD" ,"AD" ,"CTL","AD" ,"MCI","AD" ,"AD" ,"AD" ,"MCI","AD" ,"AD" ,
 "MCI","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"CTL","MCI","AD" ,"AD" ,"AD" ,"AD" ,
 "MCI","MCI","AD" ,"MCI","AD" ,"AD" ,"AD" ,"MCI","AD" ,"AD" ,"MCI","MCI",
 "CTL","AD" ,"AD" ,"MCI","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"MCI",
 "AD" ,"AD" ,"MCI","AD" ,"AD" ,"MCI","MCI","AD" ,"AD" ,"CTL","AD" ,"CTL",
 "AD" ,"AD" ,"CTL","AD" ,"CTL","AD" ,"AD" ,"MCI","MCI","MCI","AD" ,"AD" ,
 "CTL","AD" ,"AD" ,"AD" ,"MCI","AD" ,"MCI","MCI","AD" ,"AD" ,"AD" ,"AD" ,
 "AD" ,"CTL","CTL","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"MCI","AD" ,"AD" ,"AD" ,
 "AD" ,"AD" ,"MCI","AD" ,"MCI","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"MCI",
 "MCI","MCI","CTL","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"MCI","AD" ,"CTL","AD" ,
 "AD" ,"MCI","AD" ,"CTL","MCI","AD" ,"AD" ,"MCI","CTL","AD" ,"MCI","AD" ,
 "MCI","MCI","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"MCI","AD" ,"AD" ,
 "AD" ,"MCI","CTL","MCI","MCI","AD" ,"AD" ,"CTL","AD" ,"MCI","AD" ,"AD" ,
 "AD" ,"MCI","AD" ,"AD" ,"MCI","CTL","MCI","AD" ,"AD" ,"AD" ,"MCI","AD" ,
 "CTL","AD" ,"MCI","CTL","CTL","MCI","AD" ,"MCI","MCI","CTL","AD" ,"MCI",
 "AD" ,"AD" ,"AD" ,"AD" ,"MCI","AD" ,"CTL","AD" ,"MCI","CTL","MCI","AD" ,
 "AD" ,"CTL","AD" ,"MCI","MCI","MCI","AD" ,"AD" ,"MCI","AD" ,"CTL","AD" ,
 "MCI","MCI","CTL","MCI","MCI","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"MCI","AD" ,
 "CTL","AD" ,"MCI","CTL","MCI","MCI","AD" ,"CTL","AD" ,"AD" ,"AD" ,"CTL",
 "AD" ,"MCI","MCI","MCI","CTL","AD" ,"MCI","AD" ,"CTL","MCI","AD" ,"MCI",
 "AD" ,"MCI","MCI","AD" ,"CTL","AD" ,"CTL","CTL","MCI","CTL","AD" ,"AD" ,
 "AD" ,"AD" ,"MCI","AD" ,"AD" ,"AD" ,"MCI","AD" ,"MCI","AD" ,"AD" ,"MCI",
 "AD" ,"CTL","AD" ,"MCI","AD" ,"AD" ,"AD" ,"AD" ,"MCI","AD" ,"MCI","AD" ,
 "AD" ,"CTL","AD" ,"MCI","MCI","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"CTL",
 "MCI","MCI","CTL","AD" ,"AD" ,"AD" ,"MCI","AD" ,"AD" ,"CTL","AD" ,"MCI",
 "AD" ,"MCI","AD" ,"AD" ,"AD" ,"MCI","AD" ,"CTL","AD" ,"MCI","MCI","MCI",
 "MCI","AD" ,"AD" ,"MCI","AD" ,"MCI","CTL","AD" ,"AD" ,"CTL","AD" ,"AD" ,
 "MCI","MCI","AD" ,"MCI","AD" ,"AD" ,"CTL","CTL","MCI","AD" ,"MCI","AD" ,
 "AD" ,"MCI","AD" ,"MCI","AD" ,"AD" ,"AD" ,"MCI","AD" ,"AD" ,"CTL","CTL",
 "MCI","AD" ,"MCI","MCI","AD" ,"MCI","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,
 "AD" ,"CTL","AD" ,"MCI","MCI","AD" ,"MCI","AD" ,"MCI","AD" ,"AD" ,"MCI",
 "AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"CTL","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,
 "AD" ,"MCI","AD" ,"MCI","AD" ,"AD" ,"AD" ,"MCI","AD" ,"MCI","CTL","AD" ,
 "CTL","AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"AD" ,"MCI","AD" ,"AD" ,"MCI","MCI",
 "MCI","AD" ,"AD" ,"CTL","MCI" )

#######################################
# all against all linear regression
f <- factor(pData(eset)$GROUPS, levels=c("AD","MCI","CTL"))
design <- model.matrix(~0+f)
colnames(design) <- c("AD","MCI","CTL")

fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(MCI-AD,CTL-MCI,CTL-AD, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 
topTable(fit2, coef=1, adjust="BH")
# top p values are not significant
#                          logFC   AveExpr         t      P.Value  adj.P.Val
# WpHyEbnvUX1dOe9RdU -0.37672298  6.797646 -4.383855 1.452566e-05 0.05775401
# 3ep7rP6uN.E99IdO7s -0.20238061  7.460364 -3.683862 2.576398e-04 0.30231217
# NSUyuS.n3E45wKeYOo -0.17885636  9.614231 -3.578469 3.830978e-04 0.30231217

results <- decideTests(fit2)
# not hits at all...


########################################
# gene set enrichment test
# put top 50 (-1) genes into http://amp.pharm.mssm.edu/Enrichr/
# cell type works... nothing else...

# KIR2DL3 VCAN RPL12 C20orf24 NUDT1 PTPRE RNF114 LOC440055 NAPSB FKBP5 RNASE6
# CETN2 LOC645688 SLA LOC100132795 RPS14 PRMT1 CCDC12 TACC1 LOC100134648 CD1D 
# ATP6V1A EIF3D LAMP2 ARPP19 RPL10A CD79B ITGAM LPP IL18RAP STAG3L2 STK4 PSMD4 
# RPL10A PSMD4 CDK2AP1 CEBPD TM9SF2 CD302 LOC728554 TMEM188 HLA-DRA STAG3L3 
# LOC729742 NOP56 LOC728139 CRISPLD2 TP53INP1 LOC389101

########################################
# control vs AD linear regression

f <- factor(pData(eset)$GROUPS, levels=c("AD","MCI","CTL"))
design <- model.matrix(~0+f)
colnames(design) <- c("AD","MCI","CTL")

fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(CTL-AD, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 
# a = topTable(fit2, coef=1, adjust="BH",number=47231)
# write.csv(a,file='CTL_vs_AD.csv')
topTable(fit2, coef=1, adjust="BH")
#                         logFC   AveExpr         t      P.Value adj.P.Val
# KBVEl_OkmXbfQf_59c  1.0818860  7.416249  3.914995 0.0001043906 0.1037568
# cXS4CS030k.1JXMlFU  0.1679272 11.643231  3.792217 0.0001696494 0.1037568
# fNrt557tVK37urpEjc -0.1784256  7.471547 -3.761283 0.0001913390 0.1037568
# B4rSS.s4hMn11PlVHU -0.1062970 10.486157 -3.760960 0.0001915784 0.1037568
# 02tiuhr_lsCu6c8E64  0.1345540 10.521494  3.708073 0.0002348803 0.1037568
# Q0R16neBERlJ13eCYc -0.3762454  9.920894 -3.700400 0.0002418800 0.1037568
# lc17u5c4j6uS7Ne8nU -0.1857746  7.564012 -3.697792 0.0002443030 0.1037568
# 3fSd576nru6EgoOgJQ -0.3010411  9.953082 -3.684038 0.0002574673 0.1037568
# feTnq95EnpdWb6kE_c -0.2334693  8.360775 -3.669752 0.0002718438 0.1037568
# QS_tNBz.V_Hnteo6Xk -0.1872415  8.342425 -3.619981 0.0003280508 0.1037568

# no significant hits

results <- decideTests(fit2)
# not hits at all...

#######################################
# gene set enrichment at http://amp.pharm.mssm.edu/Enrichr

# top 50 probes, top 50 genes 

# One hit with GO Biological Process

# Term:
# adaptive immune response based on somatic recombination of immune receptors 
# built from immunoglobulin superfamily domains (GO:0002460)

# Overlap	P-value	Adjusted P-value	Z-score	Combined Score	Genes
# 4/57	0.00002438	0.0217	-2.162	8.28	BCL6;SLC11A1;CTSH;TLR8

# Several hits with GO Cellular Component
# Term	Overlap	P-value	Adjusted P-value	
# late endosome membrane (GO:0031902)	4/29	0.0000022011	0.0002971
# endosome membrane (GO:0010008)	6/175	0.00001137	0.000767

################################################################################
# regression on CSF measures                                                   #
################################################################################

# Load libraries
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(limma)
library(lumiHumanIDMapping)

source("code/sjnewhouse_misc_R.R")

#
load('ba3/BA3_lumi_processing_t_out/eset_final/BA3.eset_final.RData')

eset_final

# Subject_Demographics_with_chip_data_for_processing_April2015_FINAL.xlsx
# 180 subjects with CSF measures

samples_to_remove=c(
"9370786073_A","9370786073_G","9370786073_H","9370786073_L","9370786090_A",
"9370786090_B","9370786090_C","9370786090_D","9370786090_E","9370786090_G",
"9370786090_K","9370786091_A","9370786091_B","9370786091_C","9370786091_D",
"9370786091_E","9370786091_F","9370786091_G","9370786091_H","9370786091_J",
"9370786091_K","9370786091_L","9370786098_A","9370786098_B","9370786098_C",
"9370786098_D","9370786098_E","9370786098_F","9370786098_H","9370786098_K",
"9370786098_L","9371242013_A","9371242013_B","9371242013_C","9371242013_E",
"9371242013_F","9371242013_G","9371242013_J","9371242013_L","9371242063_C",
"9371242063_D","9371242063_E","9371242063_I","9371242063_J","9371242063_K",
"9371242087_B","9371242087_C","9371242087_D","9371242087_E","9371242087_F",
"9371242087_H","9371242087_I","9371242087_L","9372535005_A","9372535005_B",
"9372535005_C","9372535005_D","9372535005_E","9372535005_F","9372535005_G",
"9372535005_H","9372535005_I","9372535019_B","9372535019_C","9372535019_D",
"9372535019_E","9372535019_H","9372535019_I","9372535019_L","9402444009_A",
"9402444009_C","9402444009_D","9402444009_F","9402444009_G","9402444009_K",
"9402444009_L","9402444013_A","9402444013_B","9402444013_C","9402444013_D",
"9402444013_E","9402444013_F","9402444013_G","9402444013_H","9402444013_I",
"9402444016_B","9402444016_C","9402444016_E","9402444016_F","9402444016_G",
"9402444016_I","9402444017_A","9402444017_C","9402444017_D","9402444017_E",
"9402444017_F","9402444017_H","9402444017_I","9402444017_J","9402444017_K",
"9402444017_L","9402444020_A","9402444020_D","9402444020_F","9402444020_G",
"9402444020_H","9402444020_J","9402444020_L","9402444021_A","9402444021_F",
"9402444021_G","9402444021_J","9402444021_K","9402444021_L","9402444022_C",
"9402444022_E","9402444022_F","9402444022_H","9402444022_I","9402444022_K",
"9402444023_B","9402444023_D","9402444023_E","9402444023_F","9402444023_J",
"9402444023_K","9402444024_B","9402444024_C","9402444024_D","9402444024_F",
"9402444024_I","9402444024_J","9402444024_K","9402444024_L","9402444029_B",
"9402444029_D","9402444029_E","9402444029_F","9402444029_I","9402444029_K",
"9402444030_A","9402444030_H","9402444030_I","9402444030_L","9402444032_B",
"9402444032_C","9402444032_E","9402444032_G","9402444032_I","9402444032_J",
"9402444034_A","9402444034_B","9402444034_E","9402444034_F","9402444034_H",
"9402444034_I","9402444034_J","9402444034_L","9402444035_C","9402444035_D",
"9402444035_E","9402444035_F","9402444035_G","9402444035_I","9402444035_J",
"9402444035_K","9402444035_L","9402444040_A","9402444040_F","9402444040_H",
"9402444040_K","9402444040_L","9464921149_A","9464921149_B","9464921149_C",
"9464921149_F","9464921149_G","9464921149_H","9464921149_I","9464921149_J",
"9464921149_L","9464921165_A","9464921165_C","9464921165_D","9464921165_F",
"9464921165_G","9464921165_H","9464921165_I","9464921165_J","9464921165_K",
"9477874003_A","9477874003_E","9477874003_F","9477874003_G","9477874003_I",
"9477874003_J","9477874003_L","9477874019_B","9477874019_J","9477874019_L",
"9477874020_A","9477874020_B","9477874020_E","9477874020_F","9477874020_G",
"9477874020_J","9477874020_L","9477874026_A","9477874026_B","9477874026_D",
"9477874026_F","9477874026_G","9477874026_H","9477874026_J","9477874026_K",
"9477874026_L","9477874036_B","9477874036_D","9477874036_E","9477874036_F",
"9477874036_H","9477874036_J","9477874036_L","9477874050_A","9477874050_B",
"9477874050_C","9477874050_F","9477874050_G","9477874050_H","9477874050_J",
"9477874050_L","9477874052_D","9477874052_E","9477874052_F","9477874052_G",
"9477874052_H","9477874052_K","9477874052_L","9477874066_C","9477874066_F",
"9477874066_H","9477874066_I","9477874066_J","9477874066_K","9477874066_L",
"9477874069_B","9477874069_C","9477874069_F","9477874069_G","9477874069_H",
"9477874069_J","9477874069_K","9477874075_B","9477874075_C","9477874075_D",
"9477874075_E","9477874075_F","9477874075_G","9477874075_H","9477874079_A",
"9477874079_C","9477874079_D","9477874079_F","9477874079_G","9477874079_H",
"9477874079_J","9477874079_K","9477874079_L","9477874080_A","9477874080_B",
"9477874080_C","9477874080_D","9477874080_E","9477874080_G","9477874080_H",
"9477874080_I","9477874080_L","9477874081_A","9477874081_B","9477874081_C",
"9477874081_F","9477874081_G","9477874081_H","9477874081_J","9477874081_K",
"9477874081_L","9534190031_G","9534190031_H","9534190031_J","9534190031_L",
"9534190041_G","9534190041_H","9534190041_L","9534190085_B","9534190085_C",
"9534190085_E","9534190085_G","9534190085_H","9534190085_I","9534190110_A",
"9534190110_C","9534190110_F","9534190110_G","9534190110_H","9534190110_I",
"9534190112_A","9534190112_B","9534190112_C","9534190112_F","9534190112_G",
"9534190112_H","9534190112_I","9534190112_J","9534190112_K","9534190112_L",
"9703789027_A","9703789027_B","9703789027_C","9703789027_E","9703789027_F",
"9703789027_I","9703789027_J","9703789027_L","9703789038_A","9703789038_D",
"9703789038_E")

eset <- removeSamples_eset_lumi(eset_final,samples_to_remove)

# Subject_Demographics_with_chip_data_for_processing_April2015_FINAL.xlsx
#  columns AX and AZ 

pData(eset)$is_abnormal_abeta <- factor(c(
1,1,0,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,1,0,0,1,0,1,1,1,1,1,0,0,1,0,1,1,
0,0,0,1,1,1,0,0,0,1,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,1,1,1,0,0,1,1,0,0,0,1,
0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,1,1,0,1,1,0,0,1,0,0,0,1,1,
0,0,1,1,1,1,0,1,0,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,0,0,1,0,1,1,0,1,1,1,0,1,1,1,0,
0,1,0,1,1,0,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,1,0,0))

pData(eset)$is_abnormal_pTau <- factor(c(
0,1,1,0,1,0,1,1,0,1,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,0,0,1,
1,1,0,1,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,
0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,1,0,1,0,1,0,0,1,0,1,0,0,1,1,1,1,1,1,1,1,
1,1,0,1,0,0,1,0,0,1,0,1,1,0,1,1,1,1,0,0,1,0,0,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,
1,0,0,0,1,1,0,1,1,1,0,0,1,0,0,0,1,0,1,1,0,0,1,1))

     # 30 0 0
     # 60 0 1
     # 55 1 0
     # 35 1 1
# They don't agree that much...


#######################################
# ab42 first 
f <- factor(pData(eset)$is_abnormal_abeta)
design <- model.matrix(~f)

fit <- lmFit(eset, design)
fit <- eBayes(fit)
topTable(fit,  adjust="BH")


# nothing significant

# top 50 probes
#  [1] "r.fXqpd3Qo4r.gjzlI" "uVwv57S8_zrSxVLXLE" "Hk3WKcxF4KBL.hbWlU"
#  [4] "HRigdv4KqSbJO3j_UU" "NiFJcRLU0DE3lXTJQo" "BX7XiiR.UpMvEn1Uew"
#  [7] "9SKSKOeQiIAnkAillk" "ipRJRd5N9gdJTutBN4" "ilNVTfVyESSE3iAn0k"
# [10] "WtVe_sT_vF0.5b9I_4" "lHnn1xeQhgQ7O.Qq4k" "NotX_u7gzu5eyKuqSo"
# [13] "9jlZFNZxVjVdNCuDjE" "fdXvu97FLk4BV7RM1U" "391X5316XgEagFItAI"
# [16] "xHnn1xeQhkQ7O.Qq4k" "TpVPn6mnns.0.Helfo" "Eir7W41EKTyJOrpOiA"
# [19] "0255m6FSUH9Z5w1Lv0" "ov7lXo1Bbl4JyDNV6E" "Hm4Kedp0f9cLuvXgp4"
# [22] "BXWHeW9TUuAd9PeqwI" "NnJ1ITq_NFyKJ3pwpE" "c_gKVfo9enqhDjhehU"
# [25] "N6F6tMp4F7zoVUD3iI" "KuXKK3qj139uutHBfE" "BNUgKt5TovpLTRrI6Y"
# [28] "TK6x0JTnl18TX.6BIk" "Q._SrnlzqV3nVQ5fIg" "B7gIll1wnlQa7vRVYc"
# [31] "BqfA3QH_8osEom54mo" "rlx3rwjvf9dArlcuU4" "xfpRznNNC.WeFHi2lo"
# [34] "fuL7dKFW0eYb4R13V8" "3hfnCnAg1eASmbT3d4" "N7v5121NRa1UQviIu0"
# [37] "rrUXVUoObUUIXyTHvg" "io4r.gjzlLkgkhW56U" "NHao516.tTcef49fo0"
# [40] "xIKtHlEEkS0mrg8SuU" "Kitd7v4hDXu7lKKdbs" "W1U7e96v9G59IgiqeU"
# [43] "Bnsku19NSoLw6LfoV8" "foqeg3t9Xh2nXSt.Jc" "WqDixOtfl1U0uXfdQc"
# [46] "W6FKXGSJ4KKZ3eKA0I" "0_DHz7_krTURH3lVHk" "lqlFQi4BMfzFIokuiU"
# [49] "KcV08EPJKkelIv7Fe4" "o5RKSV6hMoJ6iQoCrc"

#######################################
# pTau
f <- factor(pData(eset)$is_abnormal_pTau )
design <- model.matrix(~f)

fit <- lmFit(eset, design)
fit <- eBayes(fit)
topTable(fit,  adjust="BH")

# one probe significant :
# QeYyFIANekn.dSk3q0 HEBP1 CHR12 13127952-13128001 adj.P.Val 0.002364033

# top 50 probes:
#  [1] "QeYyFIANekn.dSk3q0" "EfntVwoDWJpxAm.fx8" "9nRl7DVJH_oJYqo6tM"
#  [4] "fKaNGK0ojYhRjgFe3c" "0C3Sunm7rp05ew1QT0" "Eor4koo_fiVKL_LUl4"
#  [7] "0XlXNHLk57pSEjhF.4" "oA36e7sunRl7DVJ3yo" "Hg2fgVKKSI_S4lCHuU"
# [10] "oiWI0YR5HUnu4SDiFM" "rADalK5jp0ZOw1SR.o" "3ntYFxQiRAkJOTVuHk"
# [13] "oV5F615l4qunvD30tM" "i6dUKK7JK7iQzy5TuU" "Eg5j0Rf7ktaN0fppCo"
# [16] "6KO6Xk55eN_XIFJ5F4" "HeLKWh27p0Zew1SOe8" "KtnaFShHEv3c5jdSJ0"
# [19] "fggefiVKKSK_S4kiPs" "fVyA8gDQn_pSua6dOU" "QpPFWhSSKdN.ktCSug"
# [22] "ZUCXcuxfXXfHTU1S0E" "EANqUrmOnRk7DVJH_o" "ER1CPOl2pC5764iHeg"
# [25] "ivdSKCuggFOpKLKD4s" "Q.qLfhamS6dGXsN0kc" "TfOVMMNUkf6iliuJI0"
# [28] "E7rpO5ew1Scd6CeKyo" "c536E8_p_SL1fPeqPE" "itUnHeilyuog2fgVKI"
# [31] "ZklK7rq05ew3SQ.qKU" "0itKL0AUl6UE6YDXt4" "iuJ_.u1e6uhuvrj6Sg"
# [34] "l_CUovQhSXpQTogFd0" "uyj0RVF5n0fq1KXHXs" "BpQpq6Pl4nSi9CNLeg"
# [37] "6uT7nf1C9JJL3NLrlE" "lJaLrp0Rey1ScfoIeI" "lKUJ__nTlzLVJH_opQ"
# [40] "Hp05ew1ST_rrUtPGNE" "6k4vkh_cekIuCQuriw" "Q2fAXqKSK4niNm0_R0"
# [43] "NR_FUohsr7rieI2ZT4" "W4CHugDj.KGpkunRl4" "3lHUoQr0A68QSXcXtk"
# [46] "NkEC1YkOC01HXk3r60" "lkoEEjkoEp5EnSRBPU" "T3OZey1ScVKKeCSjDY"
# [49] "H57JujJ7j64tBbe6cU" "QdQghd3XurkPge9Nb0"


#######################################
# no overlapping between the two 50 probe lists

#######################################
# ab42 + pTau 
# 0 0 -> 0 
# 0 1 -> 1
# 1 0 -> 2
# 1 1 -> 3

pData(eset)$abeta_n_pTau <- factor(c(
2,3,1,2,3,2,1,1,0,1,2,3,2,3,0,1,1,2,2,2,2,1,1,3,0,0,3,1,2,2,2,2,3,1,1,2,0,2,
3,1,1,0,3,2,2,0,1,0,2,0,3,3,1,3,1,1,1,1,0,2,0,0,1,1,3,3,1,3,2,2,1,1,2,2,0,0,
0,2,0,0,0,1,1,1,1,3,3,0,2,1,2,0,1,1,1,2,1,1,2,3,0,1,0,0,3,2,1,2,2,1,1,3,1,1,
1,3,3,1,1,2,3,2,2,1,2,0,3,0,1,1,2,3,3,1,3,2,0,3,2,2,0,1,1,2,1,2,2,1,3,2,2,1,
2,3,3,0,1,2,0,2,3,1,2,3,3,1,2,2,3,2,2,0,1,2,1,1,2,2,1,1))

f <- factor(pData(eset)$abeta_n_pTau )
design <- model.matrix(~0+f)
colnames(design)
# [1] "f0" "f1" "f2" "f3"

fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(f0-f1,f1-f2,f2-f3,f0-f3,
                                  levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")

# again, nothing significant, top adj.P.Val 0.05252238

topTable(fit2, coef=1, adjust="BH", number = 50 )$ILMN_GENE

# CAT PPBP BNIP3L MAP2K2 CLC GYPE GNG11 CD36 MKNK2 RNF216L SLC20A1 C19ORF60 
# DHRS9 CORO7 GSTO1 LIPA TRABD ACAP1 HCP5 RNF167 CDK2AP1 FGL2 NCOA4 LMO2 LYZ 
# ACP1 FBXL11 PSMA1 LTA4H MIF SPG7 SKIV2L CALM2 SAMD9 NT5C3 P4HB CPNE1 CLPTM1L 
# ILF3 SAFB2 RFNG SDPR GUSBL1 FCER1A GSTO1 SLFN11 MFSD1 LOC729279 DERA HEXB

################################################################################
# gene set enrichment test 
# http://amp.pharm.mssm.edu/Enrichr/

# # GO terms, biological process
# Term	
# positive regulation of myeloid leukocyte cytokine production involved \
# in immune response (GO:0061081)	
# Overlap	P-value	Adjusted P-value	Z-score	Combined Score	Genes
# 2/11	0.00037	0.0424	-2.79333	8.8288	CD36;MIF

# # CMAP down (The Connectivity Map )
# Term	Overlap	P-value	Adjusted P-value	Z-score	Combined Score	
# dl-alpha tocopherol-1320	5/100	0.0000217	0.0164	-1.674	6.88	
# valproic acid-1155	5/100	0.0000217	0.0164	-1.6374	6.73	
# Genes
# TRABD;SAMD9;FGL2;CORO7;GNG11
# TRABD;SAMD9;FGL2;CPNE1;SKIV2L



################################################################################
# TREM2- vs R47H

# AD and MCI only 
load('BA3.eset_final.RData')
eset_final

samples_to_remove=c(
"9370786073_H","9370786090_B","9370786090_D","9370786090_E","9370786091_A",
"9370786091_C","9370786091_D","9370786091_H","9370786098_B","9370786098_C",
"9370786098_G","9370786098_I","9370786098_J","9370786098_K","9371242013_A",
"9371242013_B","9371242013_L","9371242063_C","9371242063_F","9371242063_L",
"9372535005_A","9372535019_J","9372535019_L","9402444009_D","9402444009_F",
"9402444013_A","9402444013_C","9402444013_D","9402444013_F","9402444013_H",
"9402444013_K","9402444016_E","9402444016_J","9402444016_K","9402444017_A",
"9402444017_C","9402444020_H","9402444021_E","9402444022_C","9402444022_D",
"9402444022_F","9402444022_I","9402444022_J","9402444022_K","9402444023_D",
"9402444024_C","9402444024_L","9402444029_B","9402444029_G","9402444029_H",
"9402444029_J","9402444029_K","9402444030_D","9402444030_J","9402444032_F",
"9402444032_I","9402444032_J","9402444034_B","9402444034_E","9402444035_D",
"9402444035_G","9402444035_K","9402444040_D","9402444040_J","9464921149_B",
"9464921149_C","9464921149_G","9464921165_A","9464921165_D","9464921165_H",
"9464921165_J","9464921165_L","9477874003_E","9477874003_F","9477874003_J",
"9477874003_L","9477874019_G","9477874019_L","9477874020_B","9477874020_C",
"9477874020_E","9477874020_J","9477874026_F","9477874026_K","9477874036_B",
"9477874036_D","9477874050_B","9477874050_C","9477874050_J","9477874052_B",
"9477874052_E","9477874052_L","9477874066_D","9477874066_F","9477874066_I",
"9477874066_L","9477874069_A","9477874069_B","9477874069_I","9477874075_B",
"9477874075_E","9477874075_F","9477874075_L","9477874079_B","9477874079_E",
"9477874079_F","9477874079_G","9477874079_H","9477874080_I","9477874080_J",
"9477874081_C","9477874081_D","9477874081_K","9534190031_J","9534190031_L",
"9534190041_G","9534190041_H","9534190041_K","9534190041_L","9534190085_C",
"9534190085_H","9534190085_L","9534190110_D","9534190110_F","9534190110_G",
"9534190110_H","9534190110_I","9534190112_K","9534190112_L","9703789027_C",
"9703789027_F","9703789038_A","9703789038_K")

eset=removeSamples_eset_lumi(eset_final,samples_to_remove)

dim(eset)
# Features  Samples 
#     3976      373

pData(eset)$GROUPS = c(
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","R47H"  ,"TREM2-","TREM2-",
"TREM2-","TREM2-","R47H"  ,"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","R47H"  ,"TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","R47H"  ,"TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","R47H"  ,"TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","R47H"  ,"TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","R47H"  ,"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"R47H"  ,"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"R47H"  ,"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","R47H"  ,"TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","R47H"  ,"TREM2-",
"TREM2-","TREM2-","R47H"  ,"TREM2-","TREM2-","TREM2-","TREM2-","R47H"  ,
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"R47H"  ,"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","R47H"  ,"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","TREM2-","R47H"  ,"TREM2-",
"TREM2-","TREM2-","TREM2-","TREM2-","TREM2-")

f <- factor(pData(eset)$GROUPS)
design <- model.matrix(~f)

fit <- lmFit(eset, design)
fit <- eBayes(fit)
write.table(topTable(fit,  coef="fTREM2-", adjust="BH",number=3976),
            file='TREM2-_vs_R47H_AD_n_MCI_only.txt')

################################################################################
