#!/usr/bin/Rscript
library(SNPolisher)
library(methods)

snp_file  = 'full_snp.ls'

args = commandArgs(trailingOnly=TRUE)
batch_id  = args[1]
batch_dir = args[2]
# batch_id  = 'b01'
# batch_dir = 'plates1-53/'

call_dir  = paste('/kuser/shared/data/GWAS_backup/',batch_dir,sep='')

ps_file   = paste(call_dir,'AxiomGT1.snp-posteriors.txt',sep='')
call_file = paste(call_dir,'AxiomGT1.calls.txt',         sep='')
conf_file = paste(call_dir,'AxiomGT1.confidences.txt',   sep='')
summ_file = paste(call_dir,'AxiomGT1.summary.txt',       sep='')

temp.dir  = batch_id
matr_file = paste(batch_id,'/',"metrics.txt",sep='')

refFile    = NULL
plot.prior = FALSE
priorFile  = NULL
match.cel.file.name = FALSE

################################################################################
# calculating metrics and perform classification

Ps_Metrics(
    posteriorFile      = ps_file,
    callFile           = call_file,
    output.metricsFile = matr_file
          )

Ps_Classification(
    metricsFile = matr_file ,
    output.dir  = temp.dir  ,
                 )


################################################################################
# Use a perl script to do some grab job:
# getting sub-tables for the listed SNPs only

# start reading files now

options(stringsAsFactors = FALSE)

l  <- length(summ_file)
l1 <- length(temp.dir)
l2 <- length(call_file)
if (l != l1 || l != l2)
    stop("input files need to be of the same length")
for (i in 1:l) {
    dir.create(temp.dir[i], showWarnings = F, recursive = T)
}

cmd_pidFile        <- paste("\"", snp_file,  "\"", sep = "")
cmd_summaryFile    <- paste("\"", summ_file, "\"", sep = "")
cmd_callFile       <- paste("\"", call_file, "\"", sep = "")
cmd_confidenceFile <- paste("\"", conf_file, "\"", sep = "")
cmd_posteriorFile  <- paste("\"", ps_file,   "\"", sep = "")

cmd_temp.dir       <- paste("\"", temp.dir,  "\"", sep = "")

if (is.null(refFile)) {
    cmd_refFile = "NULL"
} else {
    cmd_refFile    <- paste("\"", refFile,   "\"", sep = "")
}
if (is.null(priorFile) ||  !plot.prior ){
    cmd_priorFile = "NULL"
} else {
    cmd_priorFile  <- paste("\"", priorFile, "\"", sep = "")
}

cmd_script <-paste("\"", path.package("SNPolisher"),
                   "/Perl/visualization.pl", "\"", sep = "")

cmd <- paste(
    "perl", cmd_script,
    cmd_pidFile, cmd_summaryFile, cmd_callFile,
    cmd_confidenceFile, cmd_posteriorFile,
    cmd_temp.dir, cmd_refFile, cmd_priorFile, sep = " "
            )

sapply(cmd, system)
################################################################################