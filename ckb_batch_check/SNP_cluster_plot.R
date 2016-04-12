#!/usr/bin/Rscript
library(SNPolisher)
library(methods)

snp_file  = 'snp.ls'

# batch_id  = 'b01'
# batch_dir = 'plates1-53/'

args = commandArgs(trailingOnly=TRUE)
batch_id  = args[1]
batch_dir = args[2]

call_dir  = paste('/kuser/shared/data/GWAS_backup/',batch_dir,sep='')

ps_file   = paste(call_dir,'AxiomGT1.snp-posteriors.txt',sep='')
call_file = paste(call_dir,'AxiomGT1.calls.txt',         sep='')
conf_file = paste(call_dir,'AxiomGT1.confidences.txt',   sep='')
summ_file = paste(call_dir,'AxiomGT1.summary.txt',       sep='')

temp.dir  = batch_id
keep.temp.dir = TRUE
matr_file = paste(batch_id,'/',"metrics.txt",sep='')

geno_col=c("red","yellow","blue","gray",
           "cyan","green", "darkgreen","purple")

refFile    = NULL
sampFile   = NULL
plot.prior = FALSE
priorFile  = NULL
match.cel.file.name = FALSE

################################################################################
# calculating metrics and perform classification

#######################################
# Note we do not use the pidFile options here.
# ALL SNPs are to be classified.

to_classify_SNPs = FALSE

if (to_classify_SNPs){

    Ps_Metrics(
        posteriorFile      = ps_file,
        callFile           = call_file,
        output.metricsFile = matr_file
              )

    Ps_Classification(
        metricsFile = matr_file ,
        output.dir  = temp.dir  ,
                     )
}

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
# read the sub-tables

dd <- lapply(1:l, function(i) {
    SNPolisher:::read.snp.data(temp.dir[i], sampFile, !is.null(refFile),
                               match.cel.file.name, geno.col[6])
})

inlist <- read.delim(snp_file, header = F)[, 1]

cat("Found ", nrow(dd[[1]]$call), " of ", length(inlist),
    " requested probesets\n", sep = "")

pid <- inlist[is.element(inlist, dd[[1]]$call[, 1])]

postdata <- SNPolisher:::read.post_prior2d(
    paste(temp.dir[i],"/posterior.txt", sep = "")
                                          )

################################################################################
# start plotting now

# s = 'ps100'
for (s in pid) {
    for (i in 1:l) {
        d <- dd[[i]]
        g <- data.frame(
            called   = as.numeric(d$call[d$call[, 1] == s, -1]),
            a_signal = as.numeric(d$summary.a[d$summary.a[, 1] == s, -1]),
            b_signal = as.numeric(d$summary.b[d$summary.b[, 1] == s, -1]),
            sample = names(d$call)[-1], samps.highlight = d$samps
                       )

        if (is.null(refFile)) {
            g$reference <- rep(-1, (ncol(d$summary.a) - 1))
        } else {
            tref <- d$ref[d$ref[, 1] == s, -1]
            if (is.null(tref)) {
              g$reference <- rep(-1, (ncol(d$summary.a) - 1))
            }
            else {
              tref <- data.frame(
                         sample = names(d$ref)[-1],
                         reference = as.numeric(tref)
                                )
              g <- merge(g, tref, all.x = T)
              g$reference[is.na(g$reference)] <- -9
            }
        }

        p <- SNPolisher:::prior.for.pid(postdata, s, nclus=3)

        png(paste(temp.dir,'/',s,'.png',sep=''))
        SNPolisher:::plot.cluster(
            s, g, p,
            geno.col = geno_col, nclus = 3, type = "AvM"
                                 )
        dev.off()
    }
}

################################################################################
# the end                                                                      #
################################################################################