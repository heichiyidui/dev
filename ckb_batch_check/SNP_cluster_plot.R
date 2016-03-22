#!/usr/bin/Rscript
library(SNPolisher)
library(methods)


################################################################################
# calculating metrics etc

snp_file  = 'snp.ls'
batch_dir = 'plates1-53/'
call_dir  = paste('/kuser/shared/data/GWAS_backup/',batch_dir,sep='')

ps_file   = paste(call_dir,'AxiomGT1.snp-posteriors.txt',sep='')
call_file = paste(call_dir,'AxiomGT1.calls.txt',         sep='')
conf_file = paste(call_dir,'AxiomGT1.confidences.txt',   sep='')
summ_file = paste(call_dir,'AxiomGT1.summary.txt',       sep='')


temp.dir  = '.'
keep.temp.dir = TRUE

refFile = call_file

matr_file = "metrics.txt"

# Ps_Metrics
Ps_Metrics(
    posteriorFile      = ps_file,
    callFile           = call_file,
    output.metricsFile = matr_file
          )

# Ps_Classification
Ps_Classification(
    metricsFile = matr_file,
    output.dir  = '.',
    SpeciesType = "Diploid", GTC = FALSE
                 )

################################################################################
# start plotting now

pidFile= snp_file
summaryFile=summ_file
callFile=call_file
confidenceFile=conf_file
posteriorFile=ps_file

sampFile=NULL

plot.prior=FALSE
match.cel.file.name=FALSE
priorFile=NULL

geno.col=c("red","yellow","blue","gray",
           "cyan","green", "darkgreen","purple")


# start reading files now
options(stringsAsFactors = FALSE)

l <- length(summaryFile)
l1 <- length(temp.dir)
l2 <- length(callFile)
if (l != l1 || l != l2)
    stop("input files need to be of the same length")
for (i in 1:l) {
    dir.create(temp.dir[i], showWarnings = F, recursive = T)
}

cmd_pidFile <- paste("\"", pidFile, "\"", sep = "")
cmd_summaryFile <- paste("\"", summaryFile, "\"", sep = "")
cmd_callFile <- paste("\"", callFile, "\"", sep = "")
cmd_confidenceFile <- paste("\"", confidenceFile, "\"", sep = "")
cmd_posteriorFile <- paste("\"", posteriorFile, "\"", sep = "")
cmd_temp.dir <- paste("\"", temp.dir, "\"", sep = "")
cmd_refFile <- ifelse(is.null(refFile), "NULL", paste("\"",
    refFile, "\"", sep = ""))

cmd_priorFile <- ifelse(plot.prior & !is.null(priorFile),
    paste("\"", priorFile, "\"", sep = ""), "NULL")

cmd <- paste("perl", paste("\"", path.package("SNPolisher"),
    "/Perl/visualization.pl", "\"", sep = ""), cmd_pidFile,
    cmd_summaryFile, cmd_callFile, cmd_confidenceFile, cmd_posteriorFile,
    cmd_temp.dir, cmd_refFile, cmd_priorFile, sep = " ")
sapply(cmd, system)
dd <- lapply(1:l, function(i) {
    SNPolisher:::read.snp.data(temp.dir[i], sampFile, !is.null(refFile),
                               match.cel.file.name, geno.col[6])
})

inlist <- read.delim(pidFile, header = F)[, 1]
cat("Found ", nrow(dd[[1]]$call), " of ", length(inlist),
    " requested probesets\n", sep = "")
pid <- inlist[is.element(inlist, dd[[1]]$call[, 1])]

postdata <- SNPolisher:::read.post_prior2d(
    paste(temp.dir[i],"/posterior.txt", sep = "")
                                          )

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
        }
        else {
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

        png(paste(s,'.png',sep=''))
        SNPolisher:::plot.cluster(
            s, g, p,
            geno.col = geno.col, nclus = 3
                                 )

        dev.off()
    }
}

################################################################################
