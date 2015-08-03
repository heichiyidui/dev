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
load('ba3/BA3_lumi_processing_t_out/eset_final/BA3.eset_final.RData')

eset_final

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

eset=removeSamples_eset_lumi(eset_final,samples_to_remove)
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
# cell type works... pretty much nothing else...
KIR2DL3
VCAN
RPL12
C20orf24
NUDT1
PTPRE
RNF114
LOC440055
NAPSB
FKBP5
RNASE6
CETN2
LOC645688
SLA
LOC100132795
RPS14
PRMT1
CCDC12
TACC1
LOC100134648
CD1D
ATP6V1A
EIF3D
LAMP2
ARPP19
RPL10A
CD79B
ITGAM
LPP
IL18RAP
STAG3L2
STK4
PSMD4
RPL10A
PSMD4
CDK2AP1
CEBPD
TM9SF2
CD302
LOC728554
TMEM188
HLA-DRA
STAG3L3
LOC729742
NOP56
LOC728139
CRISPLD2
TP53INP1
LOC389101