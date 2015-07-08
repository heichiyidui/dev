################################################################################
# Batch 3 Illumina expression data analysis                                    #
################################################################################

################################################################################
# load libraries 

rm(list = ls())

library(knitr)
library(markdown)
library(shiny)

opts_chunk$set(error = TRUE, tidy = TRUE, warning = FALSE, highlight = TRUE,
    cache = TRUE, comment = NA, dev = c("png", "pdf"), fig.align = "center",
    fig.show = "asis", dpi = 92)

options(stringsAsFactors = FALSE)

# Load libraries
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(impute)
library(WGCNA)

library(gplots)
library(limma)
library(vsn)
library(MBCB)
library(lumiHumanIDMapping)
library(scatterplot3d)
library(relaimpo)
library(plyr)
library(ggplot2)
library(flashClust)

# a little housekeeping for par() and subsequent plots
def.par <- par(no.readonly = TRUE)
par(def.par)

# path to gene expression processing scripts
path_to_scripts <- "/home/kuang/dev/batch_3_expression/code"

# load 'em
source(paste(path_to_scripts, "/sjnewhouse_misc_R.R", sep = ""))

ls()

################################################################################
# set project settings and I/O

# project directory
project_dir <- "/home/kuang/dev/batch_3_expression/ba3"

# set working dir again
setwd(project_dir)

# project name
project_name <- "BA3"

# processing_date <- format(Sys.time(), "%d_%m_%Y_%s")
processing_date <- 't_out'

# output directory for lumi process and plots
out_dir <- paste(project_dir, "/", project_name, 
    "_lumi_processing_", processing_date, sep = "")

# make project pre-processing directory
make_dir_command <- paste(" if [ ! -e ", out_dir, " ]; then mkdir ", out_dir, 
    "; fi", sep = "")

system(make_dir_command)

########################################
# input files 

# genomestudio reports
in_dir <- "/home/kuang/dev/batch_3_expression/input_files/"

gs_report  <- paste(in_dir,
                     "Sample_and_Control_Probe_Profile_FinalReport.txt",sep='')
gs_probe   <- paste(in_dir,
                     "Sample_and_Control_Probe_Profile_FinalReport.txt",sep='')
gs_sample  <- paste(in_dir,
                     "sample_table_Final_Report.txt",sep='')
gs_control <- paste(in_dir,
                     "Sample_and_Control_Probe_Profile_FinalReport.txt",sep='')
anno_table <- paste(in_dir,
                     "Sample_and_Control_Probe_Profile_FinalReport.txt",sep='')

# sample information file column names MUST contain :
# Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID
pheno_file      <- paste(in_dir, "pheno_info.txt", sep="")

# batch information
tech_pheno_file <- paste(in_dir, "batch_info.txt", sep="")

# control probes 
control_table   <- paste(in_dir, 'controlprobe.txt', sep='')

########################################
# QC parameters 

# detection call rate threshold
probe_det  <- 80
sample_det <- 80

# flag for gender and sampleNetwork
sex_check <- 1  # DO THIS!! I'm not providing an option to skip this
iac_check <- 1  # DO THIS!! I'm not providing an option to skip this
iac_sd_thrs <- 2  # 2 or 3

# Model based background correction method (MLE as default) All data should
# be background correceted.  The recomended methods is MBCB (Model-based
# Background Correction for Beadarray) URL
# http://www.bioconductor.org/packages/release/bioc/html/MBCB.html
mbcb_method <- "MLE"

# Transformation method
transform_method <- "log2"  ## 'vst' # log2, vst or both

# Normalisation method
norm_method <- "rsn"  ## 'rsn' # quantile, rsn, or both

# Folks that done stuff
analyst_email <- "kuang.lin@kcl.ac.uk"
analyst_name  <- "Kuang Lin"
lab_contact_email <- "angela.k.hodges@kcl.ac.uk"
lab_contact_name  <- "Angela Hodegs"

################################################################################
# write_project_settings_to_file

# write settings to file
project_settings <- data.frame(project_dir = project_dir, 
    project_name = project_name, 
    out_dir = out_dir, gs_report = gs_report, gs_probe = gs_probe, 
    gs_sample = gs_sample, gs_control = gs_control, anno_table = anno_table, 
    pheno_file = pheno_file, tech_pheno_file = tech_pheno_file, 
    control_table = control_table,
    probe_det = probe_det, sample_det = sample_det, 
    sex_check = sex_check, iac_check = iac_check, iac_sd_thrs = iac_sd_thrs, 
    mbcb_method = mbcb_method, transform_method = transform_method, 
    norm_method = norm_method, 
    analyst_email = analyst_email, analyst_name = analyst_name, 
    lab_contact_email = lab_contact_email, lab_contact_name = lab_contact_name)

# some data wrangling
project_settings <- as.data.frame(t(project_settings))
colnames(project_settings) <- "Project_Setting"
project_settings$Project_Variable <- rownames(project_settings)
project_settings <- project_settings[, c("Project_Variable", "Project_Setting")]

# write table to out_dir
write.table(project_settings, 
    file = paste(out_dir, "/", project_name, ".project_settings.csv", sep = ""),
    row.names = FALSE, quote = FALSE, sep = ",")

# check settings
project_settings

################################################################################
# BEGIN PRE-PROCESSING Raw Expression Data                                     #
################################################################################

#######################################
# 1. read raw gene expression data

if (is.na(gs_report)) stop(" WARNING!: YOU HAVENT PROVIDED ANY DATA TO READ")

eset_raw <- lumiR(gs_report, lib.mapping = "lumiHumanIDMapping", 
    checkDupId = TRUE, detectionTh = 0.01, convertNuID = TRUE, 
    inputAnnotation = TRUE, 
    annotationColumn = c("PROBE_ID", "CHROMOSOME", "SYMBOL", "DEFINITION", 
        "ACCESSION", "ENTREZ_GENE_ID", "PROBE_TYPE", "PROBE_START", 
        "PROBE_SEQUENCE", "PROBE_CHR_ORIENTATION", "PROBE_COORDINATES",
        "CHROMOSOME", "TRANSCRIPT", "ILMN_GENE", "REFSEQ_ID", "UNIGENE_ID", 
        "SYMBOL", "PROTEIN_PRODUCT"), 
    QC = TRUE)

# use about 5G memory now
# about 20 mins

control_data=read.table(control_table, header=TRUE,sep='\t',check.names=FALSE)
head(control_data)

eset_raw <- addControlData2lumi(control_data,eset_raw)

# check it
eset_raw

# n_expression_chips
n_expression_chips <- dim(eset_raw)[2]
cat("The number of expression chips=[", dim(eset_raw)[2], "]", "\n")
# 559 subjects 
# 530 subjects after removing low detection ones
# 506 subjects after further QC

# getChipInfo
chip_id <- getChipInfo(eset_raw)$chipVersion
chip_species <- getChipInfo(eset_raw)$species
chip_probes <- getChipInfo(eset_raw)$matchedProbeNumber

chip_id       # HumanHT12_V4_0_R1_15002873_B
chip_species  # Human 
chip_probes   # 47231

#######################################
# 2. read in sample, pheno, batch and pData(eset) information

###################
# read_gs_sample_info
if (is.na(gs_sample)) 
    stop(" WARNING!: YOU HAVENT PROVIDED ANY SAMPLE INFORMATION!!!")

gs_sample

gs_sample_data <- read.table(gs_sample, 
    skip = 8, as.is = T, fill = T,  head = T, sep = "\t")

rownames(gs_sample_data) <- gs_sample_data$Sample.ID

gs_sample_data <- gs_sample_data[, names(gs_sample_data) != "X"]
gs_tech_var <- c("BIOTIN", "CY3_HYB", "HOUSEKEEPING", "LABELING", 
    "LOW_STRINGENCY_HYB",  "NEGATIVE..background.", "Noise")
sel_col <- colnames(gs_sample_data) %in% gs_tech_var

colnames(gs_sample_data) <- c(colnames(gs_sample_data[!sel_col]), 
    paste("tech.", gs_tech_var, sep = ""))

# added this as genomestudio likes to add mystery columns to the end of this
# report
n_samples <- dim(gs_sample_data)[1]  # number of rows ie samples

# save it
save(gs_sample_data, 
    file = paste(out_dir, "/", project_name, ".eset_raw.gs_sample_data.RData", 
        sep = ""))

###################
# read_pheno_info

if (is.na(pheno_file)) 
    stop(" WARNING!: YOU HAVENT PROVIDED ANY PHENOTYPE INFORMATION!!!")

pheno_file

pheno_dat <- read.table(pheno_file, as.is = T, fill = T, head = T, sep = "\t")

save(pheno_dat, 
    file = paste(out_dir, "/", project_name, ".eset_raw.pheno_dat.RData",
        sep = ""))

has_pheno_cols <- c("Sample.ID","SEX","GROUPS","TISSUE","PHENOTYPE","Study_ID"
    )  %in%  names(pheno_dat)

if ("FALSE" %in% has_pheno_cols )
     stop( paste(" WARNING!: YOU ARE MISSING ESSENTIAL SAMPLE INFORMATION!",
                 " MAKE SURE YOUR PHENO_FILE HAS: ",
                 "- Sample.ID,SEX,GROUPS,TISSUE,PHENOTYPE,Study_ID !!!"))

cat( paste("Running toupper() on PHENOTYPE, GROUP AND TISSUE variables",
           "to fix potential case issues\n"))

# fix case
pheno_dat$PHENOTYPE <- toupper(pheno_dat$PHENOTYPE)
pheno_dat$SEX <- toupper(pheno_dat$SEX)
pheno_dat$GROUPS <- toupper(pheno_dat$GROUPS)
pheno_dat$TISSUE <- toupper(pheno_dat$TISSUE)

# a quick looksee at counts
table(pheno_dat$PHENOTYPE)

table(pheno_dat$GROUPS)

table(pheno_dat$TISSUE)

table(pheno_dat$SEX)

###################
# read_batch_info

if (is.na(tech_pheno_file)) 
    stop(" WARNING!: YOU HAVENT PROVIDED ANY BATCH INFORMATION!!!")

tech_pheno_file

tech_pheno <- read.table(paste(tech_pheno_file), head = T, sep = "\t")
tech_pheno$Sentrix.Barcode <- as.character(tech_pheno$Sentrix.Barcode)
rownames(tech_pheno) <- tech_pheno$Sample.ID

# add 'tech.' to the column names
colnames(tech_pheno) <- paste("tech.", names(tech_pheno), sep = "")
colnames(tech_pheno) <- c("Sample.ID", names(tech_pheno[, -1]))

save(tech_pheno, 
    file = paste(out_dir, "/", project_name, ".eset_raw.tech_pheno.RData", 
        sep = ""))

###################
# read_pdata_info

# get pData()
eset_samples <- pData(eset_raw)

# add chip order and flad for 'has expression data' to eset_samples (pData)
eset_samples$has_expression <- 1

# making chip_oder columm
eset_samples$chip_order <- 1:dim(eset_samples)[1]
save(eset_samples, 
    file = paste(out_dir, "/", project_name, ".eset_raw.pData_samples.RData", 
        sep = ""))

###################
# quick_compare_pdata_pheno_batch_gs_sample_info

# col names
names(eset_samples)

names(gs_sample_data)

names(pheno_dat)

names(tech_pheno)

# head
head(eset_samples)

head(gs_sample_data)

head(pheno_dat)

head(tech_pheno)

# quick check these should all have the same number of rows or samples!
dim(eset_samples)
dim(gs_sample_data)
dim(pheno_dat)
dim(tech_pheno)

# Venn of Sample.ID
ex <- eset_samples$sampleID
pp <- pheno_dat$Sample.ID
tt <- tech_pheno$Sample.ID
venninput <- list(ArrayExpression = ex, Batch_Info = tt, Pheno_Info = pp)
venn(venninput)
# dev.off()

########################################
# 3. check for duplicate Study_ID

# check for duplicate Study_ID
tab_id <- table(pheno_dat$Study_ID)
tab_id_df <- as.data.frame(tab_id)
colnames(tab_id_df) <- c("Study_ID", "Freq")
dupe_samples <- subset(tab_id_df, tab_id_df$Freq >= 2)

if (dim(dupe_samples)[1] > 1) 
    cat("  WARNING!: You have duplicate Study_IDs. N=[", dim(dupe_samples)[1], 
        "]", "\n")

# show dupes
dupe_samples

# save to file
write.table(dupe_samples, 
            file = paste(out_dir, "/", project_name, ".dupe_Study_IDs.txt", 
                sep = ""), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# n_unique_study_id
n_unique_study_id <- length(unique(tab_id_df$Study_ID))
cat(" The number of unique Study_Ids=[", n_unique_study_id, "]","\n")

########################################
# 4. check Sample.IDs match in eset_samples, sample & batch info files 

cat(" Megreing pdata, pheno data, batch info and genomestudio samople data.\n")

# 4.1 merge eset_samples with pheno_dat.  Keep ALL overlaps only
eset_pheno_merge <- merge(eset_samples, pheno_dat, 
                          by.x = "sampleID", by.y = "Sample.ID")
eset_pheno_merge <- eset_pheno_merge[order(eset_pheno_merge$chip_order), ]

# check size
dim(eset_samples)
dim(eset_pheno_merge)  

# 4.2 merge eset_pheno_merge with tech_pheno
eset_pheno_batch_merge <- merge(eset_pheno_merge, tech_pheno, 
                                by.x = "sampleID", by.y = "Sample.ID")

eset_pheno_batch_merge <- 
    eset_pheno_batch_merge[order(eset_pheno_batch_merge$chip_order), ]

# check size
dim(eset_samples)
dim(eset_pheno_merge)
dim(eset_pheno_batch_merge)  

# 4.3 merge all with genomestudio final report
eset_pheno_batch_gs_merge <- merge(eset_pheno_batch_merge, gs_sample_data, 
                                  by.x = "sampleID", by.y = "Sample.ID")
eset_pheno_batch_gs_merge <- 
    eset_pheno_batch_gs_merge[order(eset_pheno_batch_gs_merge$chip_order), ]

# final look at numbers in each merged data set
dim(eset_samples)
dim(eset_pheno_merge)
dim(eset_pheno_batch_merge)
dim(eset_pheno_batch_gs_merge)

# names
names(eset_pheno_batch_gs_merge)

# looksee
head(eset_pheno_batch_gs_merge)

# quick visual check to make sure chip order is intact
# plot(eset_pheno_batch_gs_merge$chip_order, pch = 20, cex = 0.6, 
#    main = "this should be a straight line")

#######################################
# 4.4 Subset raw ExpressionSet to matched Sample.IDs & Update pData() slot.

# samples in gene expression data
samples_eset <- pData(eset_raw)$sampleID
length(samples_eset)

# samples with complete data
samples_complete_data <- eset_pheno_batch_gs_merge$sampleID
length(samples_complete_data)

# samples to remove
samples_to_remove <- (samples_eset %in% samples_complete_data) == FALSE
samples_to_remove <- pData(eset_raw)$sampleID[samples_to_remove]
length(samples_to_remove)

# rename eset_raw & save
eset_raw_preqc <- eset_raw
cat(" saving eset_raw before any qc takes place",
    "- this will be the pure un altered raw data file,", 
    "subseted to samples with pheno data.\n")

cat(" File=[", 
    paste(out_dir, "/", project_name, ".eset_raw_preqc.RData", sep = ""), 
    "]\n")

save(eset_raw_preqc,
     file = paste(out_dir, "/", project_name, ".eset_raw_preqc.RData",sep = ""))

# subset eset_raw
eset_raw <- removeSamples_eset_lumi(eset = eset_raw_preqc, 
                                    sampleRemove = samples_to_remove)
# It will warn about the subsetting didn't work on controlData. 
# Nevermind, it's done in the function.

#######################################
# 5  update pData
head(eset_pheno_batch_merge)
save(eset_pheno_batch_merge, 
     file = paste(out_dir, "/", project_name, ".eset_pheno_batch_merge.RData",
         sep = ""))

old_pdata <- pData(eset_raw)
old_pdata$old_order <- 1:dim(old_pdata)[1]

# merge with eset_pheno_batch_gs_merge
new_pdata <- merge(old_pdata, eset_pheno_batch_merge, 
                   by.x = "sampleID", by.y = "sampleID", 
                   all = TRUE, sort = FALSE)
new_pdata <- new_pdata[order(new_pdata$old_order), ]

# remove columns old_order has_expression chip_order
new_pdata <- new_pdata[, -c(2, 3, 4)]
head(new_pdata)

# update rownames
rownames(new_pdata) <- new_pdata$sampleID

# update pData slot
pData(eset_raw) <- new_pdata
dim(pData(eset_raw))

# check it
eset_raw

## n_expression_chips_with_data
n_expression_chips_with_data <- dim(eset_raw)[2]

cat(" Number of chips with complete Phenotype and Batch data=[",
    n_expression_chips_with_data,  "]",  "\n")

#######################################
# 6. Add nuID to fData,  Save updated raw ExpressionSet eset_raw

# Add nuID to fData
cat(" Add nuID to fData. \n")

fData(eset_raw)$nuID <- rownames(fData(eset_raw))
head(fData(eset_raw))

# Save updated raw ExpressionSet eset_raw
save(eset_raw, 
    file = paste(out_dir, "/", project_name, ".eset_raw.RData", sep = ""))

#######################################
# 7. Write data files to out_dir for eset_raw

# Write data files to out_dir for eset_raw
write_expression_files(eset = eset_raw, 
    outfile = paste(out_dir, "/", project_name, ".eset_raw", sep = ""))

# First time writing fData slot of eset [ ... BA3.eset_raw.fDat .txt ] 
# It crashed! 
# Removed the 4 samples with zero expression. Problem solved.

################################################################################
#  QC plots on eset_raw                                                        #
################################################################################

# basic_qc_plot_lumi

# basic_qc_plot_lumi(eset_raw)

pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_raw.basic_qc_plot_lumi.pdf", sep = ""), 
    width = 31, height = 8)
basic_qc_plot_lumi(eset_raw)
dev.off()

# outliers:
# 9371242087_K ?
# 9534190041_I ?
# 9402444035_A ?
# 9534190031_I ?

# 9477874019_C ?
# 9477874019_D ?

#####################################
# eset_samples <- pData(eset_raw)
# colnames(eset_samples) = c("sampleID","SEX","GROUPS","TISSUE"  ,"PHENOTYPE" ,
#     "Study_ID","tech.Sentrix.Barcode", "tech.SampleSection" , "tech.DOB",
#     "tech.DOC","tech.DOE","tech.RIN",  "tech.RNA_YIELD")
# pData(eset_raw) <- eset_samples

pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_raw.coloured_dendrogram_lumi.pdf", sep = ""), 
    width = 63, height = 8)
coloured_dendrogram_lumi(eset_raw)
dev.off()

# outliers:
# 9371242087_K ?
# 9534190041_I ?
# 9402444035_A ?
# 9534190031_I ?

#######################################
# pca_plot_lumi

pdf(file = paste(out_dir, "/", project_name, ".eset_raw.pca_3.pdf",
    sep = ""), width = 7, height = 7)
pca_plot_lumi(eset_raw)
dev.off()

#######################################
# need to remove the 4 subjects?  

#  samples_to_remove <- c("9371242087_K","9534190041_I",
#                         "9402444035_A","9534190031_I")
# eset_raw <- 
#     removeSamples_eset_lumi(eset = eset_raw, sampleRemove = samples_to_remove)

#######################################
# SampleNetwork Plots
pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_raw.sampleNetwork_plot_all_lumi.pdf", sep = ""), 
    width = 8, height = 8)
sampleNetwork_plot_all_lumi(eset_raw, colBy = "chip")
sampleNetwork_plot_all_lumi(eset_raw, colBy = "group")
dev.off()

# two chips to be removed 
# 9371242095 and 9464921139

# samples_to_remove <- c("9371242095_L","9371242095_F","9371242095_D",
# "9371242095_E","9371242095_H","9371242095_I","9371242095_J","9371242095_G",
# "9371242095_A","9371242095_C","9371242095_B","9371242095_K","9464921139_F",
# "9464921139_E","9464921139_I","9464921139_C","9464921139_A","9464921139_G",
# "9464921139_H","9464921139_D","9464921139_B","9464921139_J")

# eset_raw <- 
#     removeSamples_eset_lumi(eset = eset_raw, sampleRemove = samples_to_remove)

#######################################
# Save eset_raw
save(eset_raw, 
     file = paste(out_dir, "/", project_name, ".eset_raw.RData", sep = ""))
write_expression_files(eset = eset_raw, 
    outfile = paste(out_dir, "/", project_name, ".eset_raw", sep = ""))

################################################################################
# MBCB (Model-based Background Correction for Beadarray)                       #
################################################################################

#######################################
# 1. Run Model-based Background Correction for Beadarray : method = MLE  
eset_bg <- bgcor_mbcb(eset = eset_raw, 
              outfile = paste(out_dir, "/", project_name, ".eset_bg", sep = ""))

# save it
save(eset_bg, 
    file = paste(out_dir, "/", project_name, ".eset_bg.RData", sep = ""),
    compress = T)

# basic_qc_plot_lumi

pdf(file = paste(out_dir, "/", project_name, ".eset_bg.basic_qc_plot_lumi.pdf",
                 sep = ""),
    width = 61, height = 8)
basic_qc_plot_lumi(eset_bg)
dev.off()

# two  outliers  found :
# 9477874019_C 9477874019_D 
# remove them 
# samples_to_remove <-c('9477874019_C','9477874019_D')
# eset_bg <- 
#     removeSamples_eset_lumi(eset = eset_bg, sampleRemove = samples_to_remove) 

# ###################
# # pData is lost? need to re-merge them? 

# # update pData
# old_pdata <- pData(eset_bg)
# old_pdata$old_order <- 1:dim(old_pdata)[1]

# # merge with eset_pheno_batch_gs_merge
# new_pdata <- merge(old_pdata, eset_pheno_batch_gs_merge, 
#                    by.x = "sampleID", by.y = "sampleID",  sort = FALSE)
# new_pdata <- new_pdata[order(new_pdata$old_order), ]

# # remove columns old_order has_expression chip_order
# new_pdata <- new_pdata[, -c(2, 3, 4)]

# # update rownames
# rownames(new_pdata) <- new_pdata$sampleID

# # update pData slot
# pData(eset_bg) <- new_pdata
# dim(pData(eset_bg))

###################
# QC plots of eset_bg

pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg.coloured_dendrogram_lumi.pdf", sep = ""), 
    width = 61, height = 8)
coloured_dendrogram_lumi(eset_bg)
dev.off()

# pca 
pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg.pca_plot_lumi.pdf", sep = ""),
    width = 7, height = 7)
pca_plot_lumi(eset_bg)
dev.off()

# sample network 
pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg.sampleNetwork_plot_all_lumi.pdf",
    sep = ""), width = 8, height = 8)
sampleNetwork_plot_all_lumi(eset_bg, colBy = "chip")
sampleNetwork_plot_all_lumi(eset_bg, colBy = "group")
dev.off()

################################################################################
# skipped the default gender checking

# Angela did the gender checking with 4 probes: 
# ILMN_2056795  USP9Y
# ILMN_1764573  XIST
# ILMN_1755537  EIF1AY (a)
# ILMN_2228976  EIF1AY (b)
# The problematic ones are removed before anyway. 

################################################################################
# Id Detected Probes per GROUPS                                                #
################################################################################

## get neg control data from eset_bg
negativeControl <- getControlData(eset_bg)
negativeControl <- subset(negativeControl, 
                          negativeControl$controlType == "NEGATIVE")
negativeControl <- negativeControl[, c(3:dim(negativeControl)[2])]

## get neg control info mean,sd,max etc
neg_max <- apply(negativeControl, 2, max)
neg_sd <- apply(negativeControl, 2, sd)
neg_mean <- apply(negativeControl, 2, mean)
neg_2sd <- neg_mean + 2 * neg_sd

## get expression matrix
gx <- exprs(eset_bg)

## get negative bead ranges mean or max or 2SD mean of neg beads
neg_2sd <- neg_mean + 2 * neg_sd

## sweep through gx matrix to id probes greater than Mean of negative beads
## rows are probes, cols are samples

## THIS IS THE MAIN PROBE DETECTION CALCULATION
det <- sweep(gx, 1, round(neg_mean, 2), ">")

## Writing Probe Detection Calls to file #
det_tmp <- as.data.frame(det)
det_tmp$nuID <- rownames(det_tmp)
det_tmp$min_expression <- apply(gx, 1, min)
det_tmp$max_expression <- apply(gx, 1, max)
det_tmp$mean_expression <- apply(gx, 1, mean)
det_tmp$sd_expression <- apply(gx, 1, sd)

# need nuId for eset_bg
fData(eset_bg)$nuID <- rownames(fData(eset_bg))

probe_detected <- merge(fData(eset_bg), det_tmp, by.x = "nuID", by.y = "nuID",
    sort = FALSE)

cat(" Writing Probe Detection Calls to", 
    paste(out_dir, "/", project_name, ".eset_bg.probe_detected.txt", sep = ""),
    "\n")

write.table(probe_detected, 
    file = paste(out_dir, "/", project_name, ".eset_bg.probe_detected.txt",
                 sep = ""), 
    sep = "\t", row.names = FALSE, quote = FALSE)

## probe_detection counts
probe_detection <- rowSums(det)

## n samples
n_samples <- dim(gx)[2]

## probe annotations
probes_not_detected           <- probe_detection == 0
probes_detected_in_50_sample  <- probe_detection >= n_samples * 0.5
probes_detected_in_80_sample  <- probe_detection >= n_samples * 0.8
probes_detected_in_all_sample <- probe_detection == n_samples
#
probe_annotations_0_detected   <- fData(eset_bg[probes_not_detected,])
probe_annotations_50_detected  <- fData(eset_bg[probes_detected_in_50_sample,])
probe_annotations_80_detected  <- fData(eset_bg[probes_detected_in_80_sample,])
probe_annotations_100_detected <- fData(eset_bg[probes_detected_in_all_sample,])

cat(" Adding detetion call rate for all probes and samples",
    "to fData() slot for eset_bg", "\n")

fData(eset_bg)$n_detected <- probe_detection
fData(eset_bg)$n_detected_call_rate <- round(probe_detection/n_samples, 3)

fData(eset_bg)$probes_not_detected_in_any_sample <- probe_detection == 0
fData(eset_bg)$probes_detected_in_50_sample <- 
    probe_detection >= n_samples * 0.5
fData(eset_bg)$probes_detected_in_80_sample <- 
    probe_detection >= n_samples * 0.8
fData(eset_bg)$probes_detected_in_all_sample <- probe_detection == n_samples

# add min, max, mean, sd, median
fData(eset_bg)$min_expression    <- round(apply(exprs(eset_bg), 1, min), 3)
fData(eset_bg)$max_expression    <- round(apply(exprs(eset_bg), 1, max), 3)
fData(eset_bg)$mean_expression   <- round(apply(exprs(eset_bg), 1, mean), 3)
fData(eset_bg)$sd_expression     <- round(apply(exprs(eset_bg), 1, sd), 3)
fData(eset_bg)$median_expression <- round(apply(exprs(eset_bg), 1, median), 3)

## sample_detection counts
cat(" sample probe detection rate \n")

n_probes <- dim(eset_bg)[1]
sample_detection <- colSums(det)
pData(eset_bg)$n_probes_detected <- sample_detection
pData(eset_bg)$n_probes_detected_call_rate <- round(sample_detection/n_probes,3)

save(eset_bg, 
     file = paste(out_dir, "/", project_name, ".eset_bg.RData", sep = ""))

################################################################################
# good probes are of good detection rates with different groups 

## get group information from pData() slot
group_names <- unique(pData(eset_bg)$GROUPS)
group_names

groups   <- pData(eset_bg)$GROUPS
n_groups <- length(group_names)

## get expression matrix ##
gx <- exprs(eset_bg)
# get neg_mean values. Calculated previously
head(neg_mean)

################################################################################
# THIS IS THE MAIN PROBE DETECTION CALCULATION loop through each group and     #
# id probes greater than mean neg beads in X% of samples/group                 #
################################################################################

for (n in group_names) {
    cat(" Finding probes in ", probe_det/100, 
        " of sample group [", n, "]",
        " with signal intensity greather than mean intensity",
        " of the negative control beads \n ")

    group_label <- paste(n)

    sel_samples <- pData(eset_bg)$GROUPS == n

    n_samples_in_group <- dim(gx[, sel_samples])[2]

    cat(" Number of samples in group [", n, "] = ", n_samples_in_group, "\n")

    detection_matrix <- 
        sweep(gx[, sel_samples], 1, round(neg_mean, 2)[sel_samples], ">")

    group_probe_detection <- rowSums(detection_matrix) >= 
        (probe_det/100) * n_samples_in_group

    group_probe_detection_nuID <- 
        rownames(gx[group_probe_detection == TRUE, ])

    cat(" Number of probes in group [", n, "]",
        " with signal intensity greater than the mean intensity",
        " of the negative control beads = ",
        length(group_probe_detection_nuID), "\n")

    cat(" Writing probe list to ", 
        paste(out_dir, "/", project_name, ".GROUP.", group_label, 
            ".detected_probes_nuID.txt", sep = ""), 
        "\n")

    det_probes <- as.data.frame(group_probe_detection_nuID)

    colnames(det_probes) <- c("nuID")

    write.table(det_probes, 
        file = paste(out_dir, "/", project_name, ".GROUP.", group_label, 
                     ".detected_probes_nuID.txt", sep = ""), 
        row.names = FALSE, quote = FALSE, col.names = FALSE)
}

# 3680 good probes in CASE
# 3962 good probes in CTL

#################################
# writing final good probe list

cat(" writing final good probe list to", 
    paste(out_dir, "/",project_name,".detected_probes_nuID_final.txt",sep = ""),
    "\n")

system(paste("cat ", out_dir, "/", project_name, ".GROUP",
             "****.detected_probes_nuID.txt | sort | uniq > ",
             out_dir, "/", project_name, ".detected_probes_nuID_final.txt", 
             sep = ""))

# Complained 'out of memory' here. Did it manually. 
# cat BA3.GROUP****.detected_probes_nuID.txt | sort | \
#     uniq > BA3.detected_probes_nuID_final.txt
####################

good_probes <- read.table(file = 
    paste(out_dir, "/",project_name,".detected_probes_nuID_final.txt",sep = ""), 
    sep='\t',head = FALSE)
good_probes <- paste(good_probes[, 1])
n_good_probes <- length(good_probes)
cat(" Total number of good probes = ", n_good_probes, "\n")
# 3976 good probes

good_probes_annotation <- fData(eset_bg[paste(good_probes, sep = ""), ])
head(good_probes_annotation)

# add $good_probe annotation to eset_bg
fData(eset_bg)$good_probe <- fData(eset_bg)$nuID %in% good_probes

#################################
save(eset_bg, 
     file = paste(out_dir, "/", project_name, ".eset_bg.RData", sep = ""))
write_expression_files(eset = eset_bg, 
    outfile = paste(out_dir, "/", project_name,".eset_bg", sep = ""))

# Saving good probe annotations
cat(" saving good probe annotations to ", 
    paste(out_dir,"/", project_name,".detected_probes_nuID_final.***",sep = ""),
    "\n")
# or run in the shell 
# cat *GROUP*.txt | sort | uniq > BA3.detected_probes_nuID_final.txt

save(good_probes_annotation, 
    file = paste(out_dir,"/", project_name, ".detected_probes_nuID_final.RData",
                 sep = ""))

write.table(good_probes_annotation, 
    file = paste(out_dir, "/", project_name,
                 ".detected_probes_nuID_final.txt", sep = ""), 
    quote = F, sep = "\t", row.names = F)

# looksee
table(good_probes_annotation$CHROMOSOME)

####################
# plot probe detections 

pdf(file=paste(out_dir,"/",project_name,".eset_bg.detected_probe_call_rate.pdf",
               sep = ""), 
    width = 8, height = 8)
plot(good_probes_annotation$n_detected_call_rate,
     ylim = c(0, 1), pch = "*", ylab = "Call Rate",
    main = "Probe Call Rate: Detected Probes in 80% per group")
abline(h = 0.5, col = "grey", lty = 2)
abline(h = 0.8, col = "red")
dev.off()

################################################################################
# Transform and Normalise                                                      #
################################################################################
eset_bg_log2_rsn_0 <- lumiExpresso(eset_bg, 
    bg.correct = FALSE, variance.stabilize = TRUE,
    varianceStabilize.param = list(method = paste(transform_method, sep = "")),
    normalize.param = list(method = paste(norm_method, sep = "")), 
    verbose = FALSE)

# save log2 > rsn
save(eset_bg_log2_rsn_0, 
    file = paste(out_dir,"/",project_name,".eset_bg_log2_rsn_0.RData",sep = ""),
    compress = T)

################################################################################
# NOTE: might need to re-assign pData if loading from a saved R file           #
################################################################################

# write_expression_files
write_expression_files(eset = eset_bg_log2_rsn_0, 
    outfile = paste(out_dir,"/", project_name, ".eset_bg_log2_rsn_0", sep = ""))

# basic_qc_plot_lumi eset_bg_log2_rsn_0 pre-sample removal
pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg_log2_rsn_0.basic_qc_plot_lumi.pdf", sep = ""), 
    width = 11, height = 8)
basic_qc_plot_lumi(eset_bg_log2_rsn_0)
dev.off()

# coloured_dendrogram_lumi eset_bg_log2_rsn_0 pre-sample removal
pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg_log2_rsn_0.coloured_dendrogram_lumi.pdf", sep = ""),
    width = 11, height = 8)
coloured_dendrogram_lumi(eset_bg_log2_rsn_0)
dev.off()

# pca_plot_lumi eset_bg_log2_rsn_0 pre-sample removal
pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg_log2_rsn_0.pca_plot_lumi.pdf", sep = ""),
    width = 7, height = 7)
pca_plot_lumi(eset_bg_log2_rsn_0)
dev.off()

# SampleNetwork Plots eset_bg_log2_rsn_0 pre-sample removal
pdf(file = paste(out_dir, "/", project_name, 
                ".eset_bg_log2_rsn_0.sampleNetwork_plot_all_lumi.pdf",sep = ""),
    width = 8, height = 8)
sampleNetwork_plot_all_lumi(eset_bg_log2_rsn_0, colBy = "chip")
sampleNetwork_plot_all_lumi(eset_bg_log2_rsn_0, colBy = "group")
dev.off()

# age and study_id clearly matters

################################################################################
# ISA outlier removal                                                          #
################################################################################

# ISAoutliers <- basic_sampleNetworkIterate(eset = eset_bg_log2_rsn_0, 
#     col_by_chip = 0, groups = "byGroup", 
#     outfile = paste(out_dir,"/",project_name, ".eset_bg_log2_rsn_0",sep = ""),
#     IACthresh = 0.95, sd_thrs = iac_sd_thrs)

ISAoutliers <- basic_sampleNetworkIterate(eset = eset_bg_log2_rsn_0,
    col_by_chip = 0, groups = "all", 
    outfile = paste(out_dir, "/", project_name, 
                    ".eset_bg_log2_rsn_0_grp_all", sep = ""), 
    IACthresh = 0.95, sd_thrs = iac_sd_thrs)

# way too many outliers after a few rounds

# 24 outliers in the first round 

# 9370786073_D
# 9370786073_F
# 9402444013_J
# 9402444020_C
# 9402444020_I
# 9402444020_K
# 9402444021_H
# 9402444021_I
# 9402444023_A
# 9402444023_L
# 9402444024_A
# 9402444032_A
# 9402444032_L
# 9402444040_I
# 9477874019_F
# 9477874019_I
# 9477874019_K
# 9477874026_I
# 9477874036_A
# 9477874069_E
# 9534190031_K
# 9534190110_L
# 9703789027_D
# 9703789027_K

# 28 outliers in the second, don't bother

# 9370786073_J
# 9370786090_F
# 9371242087_E
# 9372535019_E
# 9402444009_I
# 9402444013_F
# 9402444013_L
# 9402444016_B
# 9402444016_F
# 9402444020_B
# 9402444020_J
# 9402444021_L
# 9402444022_J
# 9402444023_J
# 9402444024_G
# 9402444024_J
# 9402444034_G
# 9464921165_E
# 9477874020_L
# 9477874050_I
# 9477874066_C
# 9477874066_D
# 9477874075_C
# 9477874079_C
# 9534190031_I
# 9534190110_C
# 9703789027_H
# 9703789038_J

################################################################################
# PCA Batch Regressions                                                        #
################################################################################

eset_lumiN <- eset_bg_log2_rsn_0
gx <- exprs(eset_lumiN)
gx <- t(gx)

pca_gx <- prcomp(gx)$x
pca_gx <- pca_gx[, "PC1"]

# PHENOTYPES FOR REGRESSIONS
cat(" setting up phenotypes PC1,PHENOTYPE & GROUPS for regressions \n")

PC1 <- as.numeric(pca_gx)
PHENOTYPE <- as.numeric(as.factor(toupper(pData(eset_lumiN)$PHENOTYPE))) 
# toupper() called because of pesky 'case' issues
GROUPS <- as.numeric(as.factor(toupper(pData(eset_lumiN)$GROUPS)))  
# toupper() called because of pesky 'case' issues

# df <- cbind(tech_batch,PC1,PHENOTYPE,GROUPS) 
# df_z <- apply(df,2,as.factor)
# df_z <- apply(df_z,2,as.numeric)

pdata <- pData(eset_lumiN)
pdata <- as.data.frame(pdata[,
    c("tech.Sentrix.Barcode","tech.SampleSection","SEX","Study_ID",
      "tech.AGE","tech.DOE","tech.RIN","tech.RNA_YIELD")])

# LK: got some 'NA' in AGE, RIN and RNA_YIELD
# replace them with the mean
pdata$tech.AGE[6  ]=77
pdata$tech.AGE[15 ]=77
pdata$tech.AGE[23 ]=77
pdata$tech.AGE[30 ]=77
pdata$tech.AGE[37 ]=77
pdata$tech.AGE[75 ]=77
pdata$tech.AGE[82 ]=77
pdata$tech.AGE[110]=77
pdata$tech.AGE[117]=77
pdata$tech.AGE[132]=77
pdata$tech.AGE[140]=77
pdata$tech.AGE[141]=77
pdata$tech.AGE[202]=77
pdata$tech.AGE[213]=77
pdata$tech.AGE[234]=77
pdata$tech.AGE[240]=77
pdata$tech.AGE[271]=77
pdata$tech.AGE[284]=77
pdata$tech.AGE[297]=77
pdata$tech.AGE[298]=77
pdata$tech.AGE[313]=77
pdata$tech.AGE[316]=77
pdata$tech.AGE[321]=77
pdata$tech.AGE[323]=77
pdata$tech.AGE[348]=77
pdata$tech.AGE[355]=77
pdata$tech.AGE[375]=77
pdata$tech.AGE[378]=77
pdata$tech.AGE[398]=77
pdata$tech.AGE[407]=77
pdata$tech.AGE[431]=77
pdata$tech.AGE[439]=77
pdata$tech.AGE[468]=77
pdata$tech.AGE[495]=77

pdata$tech.RIN[75]=8.4
pdata$tech.RIN[323]=8.4
pdata$tech.RIN[407]=8.4
pdata$tech.RNA_YIELD[443] = 834.95
# random_shuffle_line to select some random dates
pdata$tech.DOE[75]='28/11/13'
pdata$tech.DOE[323]='28/11/13'
pdata$tech.DOE[407]='09/01/14'

# used as factors
pdata$SEX                  <- as.factor(pdata$SEX)
pdata$Study_ID             <- as.factor(pdata$Study_ID)
pdata$tech.Sentrix.Barcode <- as.factor(pdata$tech.Sentrix.Barcode)
pdata$tech.SampleSection   <- as.factor(pdata$tech.SampleSection)

pdata$tech.DOE             <- as.factor(pdata$tech.DOE )

tech_batch <- pdata

batch_var_names <- names(tech_batch)
batch_var_names

multivariate_model_terms <- paste(batch_var_names, collapse = "+")

# multivariate_model_terms <- paste(
#     "tech.Sentrix.Barcode+tech.SampleSection",
#     "tech.RIN+tech.RNA_YIELD+tech.DOE+SEX",sep='+')

########################
# PC1 is run last START LINEAR REGRESSION 
########################
# for(pheno in c('PHENOTYPE','PC1') ) {
for (pheno in c("PC1")) {
    # pheno <- 'PC1'
    
    multivariate_model <- paste(pheno, multivariate_model_terms, sep = "~")
    multivariate_model <- as.formula(multivariate_model)
    
    # multivariate lm
    cat(" running full multivariate models ", pheno, " ~ multivariate_model\n")
    lm_batch <- lm(multivariate_model, data = tech_batch,na.action=na.fail)
    summary(lm_batch)
    # RSQUARED summary lm
    lm_r2 <- round(summary(lm_batch)$adj.r.squared, 3)
    lm_r2 
    
    # summary lm
    summary_lm_batch <- summary(lm_batch)$coef
    summary_lm_batch <- as.data.frame(summary_lm_batch)
    summary_lm_batch$terms <- rownames(summary_lm_batch)
    summary_lm_batch$significant <- 
        ifelse(summary_lm_batch$"Pr(>|t|)" <= 0.05, 1, 0)
    summary_lm_batch$model_rsq <- lm_r2

    # save summary lm
    write.table(summary_lm_batch, 
        file = paste(out_dir, "/", project_name, ".eset_bg_log2_rsn.",
                    pheno, ".multivariate_model_batch_variables.csv", sep = ""),
        row.names = FALSE,
        quote = FALSE, sep = ",")

    # multivariate ANOVA
    anova_lm_batch <- anova(lm_batch)
    anova_lm_data <- as.data.frame(anova_lm_batch)
    anova_lm_data$terms <- rownames(anova_lm_data)
    anova_lm_data <- subset(anova_lm_data, anova_lm_data$terms != "Residuals")

    ## plot ANOVA P
    par(mar = c(15, 5, 4, 2))
    barplot(-log10(anova_lm_data$"Pr(>F)"), srt = 45, las = 3, 
        names = c(anova_lm_data$terms),
        ylab = "ANOVA -log10(P)", 
        main = paste(pheno, "~multivariate_model. R2=",lm_r2, sep = ""), 
        cex.names = 0.8, cex.main = 0.8, cex.lab = 1)
    abline(h = -log10(0.05), col = "blue")
    abline(h = -log10(0.05/dim(anova_lm_data)[1]), col = "red")

    ## plot to pdf
    pdf(file = paste(out_dir, "/", project_name, ".eset_bg_log2_rsn.", pheno,
                     ".ANOVA_multivariate_model_batch_variables.pdf", sep = ""),
        width = 8, height = 6)
    
    par(mar = c(15, 5, 4, 2))
    barplot(-log10(anova_lm_data$"Pr(>F)"), 
        srt = 45, las = 3, names = c(anova_lm_data$terms),
        ylab = "ANOVA -log10(P)", 
        main = paste(pheno, "~ multivariate_model. R2=", lm_r2, sep = ""), 
        cex.names = 0.8, cex.main = 0.8, cex.lab = 1)
    abline(h = -log10(0.05), col = "blue")
    abline(h = -log10(0.05/dim(anova_lm_data)[1]), col = "red")
    dev.off()

    ## save ANOVA
    write.table(anova_lm_data, 
        file = paste(out_dir, "/", project_name, ".eset_bg_log2_rsn.", pheno, 
                     ".ANOVA_multivariate_model_batch_variables.csv", sep = ""),
        row.names = FALSE, quote = FALSE, sep = ",")

    # are there any sig ANOVA terms
    min_anova_p <- min(anova_lm_data$"Pr(>F)")

    # if sig ANOVA terms then do this:-
    if (min_anova_p <= 0.05) {
        sig_anova_lm_data <- subset(anova_lm_data, anova_lm_data$"Pr(>F)" <=
            0.05)
        most_sig_term <- sig_anova_lm_data$terms[sig_anova_lm_data$"Pr(>F)" ==
            min_anova_p]
        sig_terms <- sig_anova_lm_data$terms
        sig_terms <- paste(sig_terms, collapse = " + ")
        cat(" WARNING!: SIGNIFICANT ASSOCIATION BETWEEN [ ", pheno, " ] ~ [",
            sig_terms, " ]. R2=", 
            lm_r2, " [", most_sig_term, "] MIN_P=", min_anova_p,
            ". YOU MAY WANT TO CORRECT FOR THIS BEFORE THE FINAL ANLYSIS \n ")
        
        ######################## STEP find best terms#
        cat(" Finding a best model by AIC in a Stepwise Algorithm \n")
        step_lm_batch <- stepAIC(lm_batch, direction = "both")
        # summary step lm
        summary_step_lm_batch <- summary(step_lm_batch)$coef
        summary_step_lm_batch <- as.data.frame(summary_step_lm_batch)
        summary_step_lm_batch$terms <- rownames(summary_step_lm_batch)
        # save summary step lm
        write.table(summary_step_lm_batch, 
            file = paste(out_dir, "/", project_name,
             ".eset_bg_log2_rsn.stepAIC_multivariate_model_batch_variables.csv",
             sep = ""), 
            row.names = FALSE, quote = FALSE, sep = ",")

        ## anova step #
        anova_step_lm_batch <- anova(step_lm_batch)
        anova_data <- as.data.frame(anova_step_lm_batch)
        anova_data$terms <- rownames(anova_data)
        anova_data <- subset(anova_data, anova_data$terms != "Residuals")
        # best model
        best_model <- paste(anova_data$terms, collapse = " + ")
        best_model <- as.formula(paste(pheno, "~", best_model, sep = ""))
        best_model

        # RSQ anova step
        anova_step_r2 <- round(summary(step_lm_batch)$adj.r.squared, 3)

        ## plot
        par(mar = c(15, 5, 4, 2))
        barplot(-log10(anova_data$"Pr(>F)"), las = 3, 
            names = c(anova_data$terms),
            ylab = "ANOVA -log10(P)", 
            main = paste(step_lm_batch$call[2], " R2=",anova_step_r2, sep = ""),
            cex.names = 0.8, cex.main = 0.6, cex.lab = 1)
        abline(h = -log10(0.05), col = "blue")
        abline(h = -log10(0.05/dim(anova_data)[1]), col = "red")

        ## plot ANOVA step P #
        pdf(file = paste(out_dir, "/", project_name, ".eset_bg_log2_rsn.",pheno,
                ".stepANOVA_multivariate_model_batch_variables.pdf", sep = ""),
            width = 8, height = 6)
        par(mar = c(15, 5, 4, 2))
        barplot(-log10(anova_data$"Pr(>F)"), las = 3, 
            names = c(anova_data$terms),
            ylab = "ANOVA -log10(P)", 
            main = paste(step_lm_batch$call[2], " R2=",anova_step_r2, sep = ""),
            cex.names = 0.8, cex.main = 0.6, cex.lab = 1)
        abline(h = -log10(0.05), col = "blue")
        abline(h = -log10(0.05/dim(anova_data)[1]), col = "red")
        dev.off()

        ## save stepANOVA
        write.table(anova_data, 
            file = paste(out_dir, "/", project_name, ".eset_bg_log2_rsn.",pheno,
                ".stepANOVA_multivariate_model_batch_variables.csv", sep = ""),
            row.names = FALSE, quote = FALSE, sep = ",")
        cat(" BEST MODEL BASED ON stepAIC [", pheno, "] ~ [", 
            paste(best_model)[3], "] ", "\r", "\n")
    } else {
        cat(" NICE!: NO SIGNIFICANT ASSOCIATION BETWEEN [ ", pheno, 
            " ] ~ [ BATCH PHENOTYPES ] \n")
    }
}
#  WARNING!: SIGNIFICANT ASSOCIATION BETWEEN 
# [  PC1  ] ~ [ tech.Sentrix.Barcode + tech.SampleSection + Study_ID 
#             + tech.DOE + tech.RIN  ]. 
# R2= 0.328  [ Study_ID ]
# MIN_P= 1.445193e-20 . 
# YOU MAY WANT TO CORRECT FOR THIS BEFORE THE FINAL ANLYSIS 
#   Finding a best model by AIC in a Stepwise Algorithm 
# Start:  AIC=3378.1
# PC1 ~ tech.Sentrix.Barcode + tech.SampleSection + SEX + Study_ID + 
#     tech.AGE + tech.DOE + tech.RIN + tech.RNA_YIELD

#                         Df Sum of Sq    RSS    AIC
# - tech.DOE             150    114361 285945 3336.5
# - tech.Sentrix.Barcode  45     23920 195505 3354.1
# - tech.AGE               1        66 171651 3376.3
# - tech.RNA_YIELD         1        96 171681 3376.4
# <none>                               171584 3378.1
# - SEX                    1       706 172290 3378.2
# - tech.SampleSection    11     11533 183117 3389.0
# - Study_ID               4     13636 185221 3408.8
# - tech.RIN               1     17622 189206 3425.6

# Step:  AIC=3301.46
# PC1 ~ tech.SampleSection + Study_ID + tech.RIN + tech.RNA_YIELD
#                         Df Sum of Sq    RSS    AIC
# <none>                               321265 3301.5
# - tech.RNA_YIELD         1      1700 322965 3302.1
# + tech.AGE               1       663 320601 3302.4
# + SEX                    1       311 320954 3303.0
# - tech.SampleSection    11     17547 338812 3306.4
# + tech.Sentrix.Barcode  45     35160 286105 3332.8
# - tech.RIN               1     26452 347717 3339.5
# + tech.DOE             150    125301 195964 3351.3
# - Study_ID               4     40964 362229 3354.2

#  BEST MODEL BASED ON stepAIC 
# [ PC1 ] ~ [ tech.SampleSection + Study_ID + tech.RIN + tech.RNA_YIELD ]  

################################################################################
# Batch Correction using linear models                                         #
################################################################################

if (pheno != "PC1") 
    stop(" WARNING!: model terms are not from the PC1 ~ batch regressions")

# get gene expressuion
gx <- exprs(eset_bg_log2_rsn_0)
n_probes <- dim(gx)[1]
gx <- t(gx)
dim(gx)
# 506 47231

# probe nuID names
new_probe_names <- paste("p_", colnames(gx), sep = "") 
 # ADD p as some nuID start with a number
# probe_names <- colnames(gx)
head(new_probe_names)
colnames(gx) <- new_probe_names
gx[1:10, 1:10]

adj_gx <- gx * 0
adj_gx[1:10, 1:10]

# get batch phenos
batch_pheno <- tech_batch[, anova_data$terms]
batch_pheno <- data.frame(cbind(batch_pheno, gx))
# this is the data for the regression

eset_bg_log2_rsn_regression_input <- batch_pheno

save(eset_bg_log2_rsn_regression_input, 
    file = paste(out_dir, "/", project_name,
                 ".eset_bg_log2_rsn.regression_input.RData", sep = ""))

########################################
# loop through each probe and adjust for sig batches
pn <- 1
#
for (probe in new_probe_names) {
    # probe <- "p_Ku8QhfS0n_hIOABXuE"

    lm_model <- as.formula(paste(probe, "~", best_model[3], sep = ""))
    lm_probe <- lm(lm_model, data = eset_bg_log2_rsn_regression_input)
    rsq <- round(summary(lm_probe)$adj.r.squared, 3)

    residual_probe <- lm_probe$residual
    mean_probe_level <- mean(batch_pheno[, probe])
    adjusted_probe_level <- residual_probe + mean_probe_level

    adj_gx[, probe] <- adjusted_probe_level

    ### cat(' Progress: ',pn,' : ',round(pn/n_probes,3),'\r')
    sink(file = paste(out_dir, "/", project_name, 
                      ".eset_bg_log2_rsn.lm_probe_progress_rsq.txt", sep = ""),
         append = TRUE)
    cat(" doing [", probe, "] ~ [", paste(best_model[3]), "].RSQ=",
        rsq, ". Progress:", round(pn/n_probes, 3), "\n")
    sink()
    pn <- pn + 1
}
# This loop took about 10 hours!

# update names and transform back to probe x sample matrix
adj_gx <- t(adj_gx)
rownames(adj_gx) <- rownames(exprs(eset_bg_log2_rsn_0))
adj_gx[1:10, 1:10]

eset_bg_log2_rsn_adj_gx <- adj_gx
save(eset_bg_log2_rsn_adj_gx, 
    file = paste(out_dir, "/", project_name, ".eset_bg_log2_rsn_adj_gx.RData",
                 sep = ""))

# make new eset and replace exprs() matrix with new batch adjusted data
eset_bg_log2_rsn_adj <- eset_bg_log2_rsn_0
exprs(eset_bg_log2_rsn_adj) <- adj_gx

# save eset_bg_log2_rsn_adj
save(eset_bg_log2_rsn_adj, 
     file = paste(out_dir, "/", project_name, ".eset_bg_log2_rsn_adj.RData",
                  sep = ""),
     compress = T)
# write_expression_files eset_bg_log2_rsn_adj
write_expression_files(eset = eset_bg_log2_rsn_adj, 
    outfile = paste(out_dir,"/",project_name,".eset_bg_log2_rsn_adj", sep = ""))

# QC plots of the adjusted eset 
pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg_log2_rsn_adj.basic_qc_plot_lumi.pdf",sep = ""), 
    width = 31, height = 8)
basic_qc_plot_lumi(eset_bg_log2_rsn_adj)
dev.off()

pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg_log2_rsn_adj.coloured_dendrogram_lumi.pdf",sep = ""),
    width = 31, height = 8)
coloured_dendrogram_lumi(eset_bg_log2_rsn_adj)
dev.off()

pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg_log2_rsn_adj.pca_plot_lumi.pdf",sep = ""), 
    width = 7, height = 7)
pca_plot_lumi(eset_bg_log2_rsn_adj)
dev.off()

pdf(file = paste(out_dir, "/", project_name, 
                 ".eset_bg_log2_rsn_adj.sampleNetwork_plot_all_lumi.pdf",
                 sep = ""), 
    width = 8, height = 8)
sampleNetwork_plot_all_lumi(eset_bg_log2_rsn_adj, colBy = "chip")
sampleNetwork_plot_all_lumi(eset_bg_log2_rsn_adj, colBy = "group")
dev.off()

################################################################################
# Create Final QC'd Expression data set                                        #
################################################################################
# subset to good probes 
eset_final <- eset_bg_log2_rsn_adj[good_probes,]
eset_final


save(eset_final, file = paste(out_dir, "/", project_name, 
                              ".eset_final.RData", sep = ""), 
    compress = T)

write_expression_files(eset = eset_final, 
    outfile = paste(out_dir, "/", project_name,".eset_final", sep = ""))

################################################################################
# Clean up                                                                     #
################################################################################

system(paste(" mkdir ", out_dir, "/eset_raw", sep = ""))
system(paste(" mkdir ", out_dir, "/eset_bg", sep = ""))
system(paste(" mkdir ", out_dir, "/eset_bg_log2_rsn", sep = ""))
system(paste(" mkdir ", out_dir, "/eset_final", sep = ""))
system(paste(" mkdir ", out_dir, "/XIST_Gender_checks", sep = ""))
system(paste(" mkdir ", out_dir, "/detected_probes", sep = ""))
system(paste(" mkdir ", out_dir, "/SampleNetwork", sep = ""))
system(paste(" mkdir ", out_dir, "/batch_regressions/", sep = ""))

# XIST
system(paste(" mv -v ", out_dir, "/", "*.XIST* ",
             out_dir,"/XIST_Gender_checks/", sep = ""))

# detected probes
system(paste(" mv -v ", out_dir, "/", "*.detected_probes* ",
             out_dir, "/detected_probes/", sep = ""))

# SampleNetwork
system(paste("  mv -v ", out_dir, "/", "*ampleNetwork* ", 
             out_dir, "/SampleNetwork/", sep = ""))

# batch regressions multivariate
system(paste(" mv -v ", out_dir, "/", "*multivariate* ", 
             out_dir, "/batch_regressions/", sep = ""))

# eset_****
system(paste(" mv -v ", out_dir, "/", "*.eset_raw.* ", 
             out_dir, "/eset_raw/",  sep = ""))
system(paste(" mv -v ", out_dir, "/", "*.eset_bg.* ", 
             out_dir, "/eset_bg/", sep = ""))

system(paste(" mv -v ", out_dir, "/", "*.eset_bg_log2_rsn* ", 
             out_dir, "/eset_bg_log2_rsn/", sep = ""))
system(paste(" mv -v ", out_dir, "/", "*.eset_final.* ", 
             out_dir, "/eset_final/", sep = ""))

################################################################################
# Create and Write Project Summary                                             #
################################################################################

n_gender_fails=0

project_summary <- data.frame(
    project_dir = project_dir, 
    project_name = project_name,
    out_dir = out_dir, gs_report = gs_report, gs_probe = gs_probe, 
    gs_sample = gs_sample, gs_control = gs_control, 
    anno_table = anno_table, pheno_file = pheno_file, 
    tech_pheno_file = tech_pheno_file, probe_det = probe_det, 
    sample_det = sample_det, 
    sex_check = sex_check, iac_check = iac_check, iac_sd_thrs = iac_sd_thrs,
    mbcb_method = mbcb_method, transform_method = transform_method, 
    norm_method = norm_method,
    analyst_email = analyst_email, analyst_name = analyst_name, 
    lab_contact_email = lab_contact_email, lab_contact_name = lab_contact_name,
    chip_id = chip_id, chip_species = chip_species,  chip_probes = chip_probes,
    n_expression_chips = n_expression_chips, 
    n_unique_study_id = n_unique_study_id,
    n_expression_chips_with_data = n_expression_chips_with_data, 
    n_gender_fails = n_gender_fails,
    n_unique_study_id_gender_fails = n_unique_study_id_gender_fails, 
    gender_concordance = gender_concordance,
    n_good_probes = n_good_probes, n_outliers = length(outlier_samples), 
    n_unique_study_id_outliers = n_unique_study_id_outliers,
    n_good_samples = dim(eset_bg_log2_rsn_0)[2], 
    n_unique_study_id_good_samples = 
        length(unique(pData(eset_bg_log2_rsn_0)$Study_ID)))

# some data wrangling
project_summary <- as.data.frame(t(project_summary))
colnames(project_summary) <- "Project_Setting"
project_summary$Project_Variable <- rownames(project_summary)
project_summary <- project_summary[, c("Project_Variable", "Project_Setting")]

# write table to out_dir
write.table(project_summary, 
    file = paste(out_dir, "/", project_name, ".project_summary.csv", sep = ""), 
    row.names = FALSE, quote = FALSE, sep = ",")

# looksee
project_summary

################################################################################
# THE END                                                                      #
################################################################################

