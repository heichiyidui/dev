source('../batch_3_expression/code/sjnewhouse_misc_R.R')

# read the protein activity data 
exprs   <- read.csv("t_exprs.csv", header = TRUE)

# histogram 
pdf('hist.pdf',width=12)
hist(exprs)
dev.off()

# run PCA
pca_gx  <- prcomp(t(exprs))
pca_1_2 <- pca_gx$rotation[1:111,1:2]

# read the demographics data
pheno  <- read.csv('t_pheno.csv',header = TRUE)

# plot the pheno types against PC1 and PC2
pdf('pca.pdf')
age_color    <- labels2colors(pheno$Age_at_sampling__of_TREM2._to_match>77)
trem2_color  <- labels2colors(pheno$is_TREM2)
plate_color  <- labels2colors(pheno$Plate_location_for_assay_Duplicate_2)
diagno_color <- labels2colors(pheno$Diagnosis_at_sampling_TREM2)
gender_color <- labels2colors(pheno$Gender_of_TREM2)
apoe_color   <- labels2colors(pheno$APOE_genotype)

plot(pca_gx$rotation[,1:2],pch=21,main=('AGE'),      bg=age_color)
plot(pca_gx$rotation[,1:2],pch=21,main=('is_TREM2+'),bg=trem2_color)
plot(pca_gx$rotation[,1:2],pch=21,main=('plate'),    bg=plate_color)
plot(pca_gx$rotation[,1:2],pch=21,main=('Diagnosis'),bg=diagno_color)
plot(pca_gx$rotation[,1:2],pch=21,main=('Gender'),   bg=gender_color)
plot(pca_gx$rotation[,1:2],pch=21,main=('Apoe'),     bg=apoe_color)

dev.off()

########################################
# T test trem2-/+  
res=t.test(exprs$FGF.basic~pheno$is_TREM2.); res$p.value 
# ...
# probe_id     t_test_p(trem2-/+)     
# EGF          0.1214684
# IFN.a        0.1307003
# IL.1b        0.1687089
# MCP.1        0.2067498
# IL.2         0.2602391
# IL.7         0.2986386
# VEGF         0.310784
# IFN.y        0.319765
# IL.13        0.3357263
# G.CSF        0.3468701
# GM.CSF       0.3488988
# IL.5         0.4023042
# Eotaxin      0.4122602
# IL.6         0.4173261
# IL.17        0.4468384
# MIP.1a       0.4532591
# IL.8         0.4661743
# TNF.a        0.5284058
# IL.4         0.5467611
# FGF.basic    0.586326
# IL.10        0.5918514
# IL.1RA       0.6066359
# MIP.1b       0.6638101
# IL.2R        0.8423331
# IL.12        0.863311
# IL.15        0.8646292
# MIG          0.8683986
# IP.10        0.8962887
# Rantes       0.9551692
# HGF          0.9831364


########################################
# T test AD/Other

is_AD <- pheno$Diagnosis_at_sampling_TREM2. == "AD"
res=t.test(exprs$FGF.basic ~ is_AD); res$p.value 
# ...

# IL.1RA    0.04046155
# Rantes    0.05966398
# EGF       0.09256972
# MIP.1b    0.1080252
# Eotaxin   0.1187788
# TNF.a     0.1352726
# IL.15     0.1608439
# IL.4      0.1786776
# IL.2R     0.1945885
# IL.7      0.2263806
# IL.8      0.2473178
# HGF       0.2538089
# G.CSF     0.2699736
# IFN.y     0.3354629
# IL.10     0.4208607
# IP.10     0.4704991
# IL.2      0.5137977
# MIP.1a    0.5296729
# FGF.basic 0.5456953
# IFN.a     0.5467529
# IL.5      0.6455544
# IL.17     0.6572386
# IL.12     0.6872026
# IL.13     0.7121233
# GM.CSF    0.768296
# IL.1b     0.7743325
# MCP.1     0.7813397
# MIG       0.8288261
# VEGF      0.93412
# IL.6      0.9908608


########################################
# T test TREM2- vs 140G>A_R47H_Hardy

pheno$TREM2._variant
f1 <-  factor(pheno$TREM2._variant,levels=c('TREM2-','140G>A_R47H_Hardy'))

res=t.test(exprs$FGF.basic ~ f1 ) ; res$p.value 
# ... 
# probe     p.value 
# IFN.y     0.0214368
# IL.17     0.04438396
# G.CSF     0.04599267
# VEGF      0.06093905
# MCP.1     0.07134537
# IFN.a     0.07585516
# IL.8      0.08110076
# IL.7      0.09054079
# EGF       0.09628497
# IL.5      0.1255056
# GM.CSF    0.1329183
# MIP.1a    0.1342651
# TNF.a     0.1524458
# IL.4      0.1718852
# IL.6      0.1830897
# IL.2R     0.2394553
# Eotaxin   0.2526601
# FGF.basic 0.2540692
# IL.13     0.2630358
# MIP.1b    0.295966
# IL.15     0.3019291
# IL.10     0.3084335
# IL.1b     0.3298444
# IP.10     0.3403181
# IL.1RA    0.3793772
# MIG       0.3945722
# HGF       0.4066752
# IL.12     0.4329994
# IL.2      0.4438556
# Rantes    0.5280287

# Slightly better. Still not significant after multiple teset correction.
################################################################################