#!/usr/bin/Rscript
library(SNPolisher)
library(methods)

################################################################################
# calculating metrics etc

# Ps_Metrics
Ps_Metrics(
    posteriorFile="example_1_buffalo/AxiomGT1.snp-posteriors.txt",
    callFile="example_1_buffalo/AxiomGT1.calls.txt",
    output.metricsFile="example_1_buffalo/Output/metrics.txt"
)

# Ps_Classification
Ps_Classification(
    metricsFile="example_1_buffalo/Output/metrics.txt",
    ps2snpFile="example_1_buffalo/ps2snp.txt",
    output.dir="example_1_buffalo/Output",
    SpeciesType="Diploid", GTC=FALSE
)

ps.performance <- read.table(
    "example_1_buffalo/Output/Ps.performance.txt",
    header=T
)

ps.performance[1:5,]
ps.performance[1:5,c(1:2,16:17)]

# pulling out SNPs for PolyHighRes and MonoHighRes
ps.performance$ConversionType[1:5]
ps.performance[1:5,16]

names(ps.performance)
unique(ps.performance$ConversionType)

converted <- ps.performance[
    ps.performance$ConversionType %in%
        c("PolyHighResolution","MonoHighResolution"),1
]

write.table(
    converted,
    file="example_1_buffalo/Output/converted.txt",
    sep="\t",
    row.names=FALSE,
    col.names="probeset_id",
    quote=FALSE
)

################################################################################
# start plotting now

pidFile="example_1_buffalo/Output/PolyHighResolution.ps"
summaryFile="example_1_buffalo/AxiomGT1.summary.txt"
callFile="example_1_buffalo/AxiomGT1.calls.txt"
confidenceFile="example_1_buffalo/AxiomGT1.confidences.txt"
posteriorFile="example_1_buffalo/AxiomGT1.snp-posteriors.txt"
sampFile="example_1_buffalo/sample_list.txt"
temp.dir="example_1_buffalo/Output/Temp"
keep.temp.dir=FALSE
refFile="example_1_buffalo/ref.txt"
plot.prior=FALSE
match.cel.file.name=FALSE
priorFile=NULL

geno.col=c("red","dodgerblue2","blue","gray",
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

inlist <- read.delim(pidFile, header = T)[, 1]
cat("Found ", nrow(dd[[1]]$call), " of ", length(inlist),
    " requested probesets\n", sep = "")
pid <- inlist[is.element(inlist, dd[[1]]$call[, 1])]

postdata <- SNPolisher:::read.post_prior2d(
    paste(temp.dir[i],"/posterior.txt", sep = "")
                                          )

nclus = 3
# s = 'ps100'
for (s in pid[1:10]) {
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

        p <- SNPolisher:::prior.for.pid(postdata, s, nclus)

        png(paste(s,'.png',sep=''))
        SNPolisher:::plot.cluster(
            s, g, p, "called", "AvM",
            geno.col = geno.col, nclus = nclus
                                 )

        dev.off()
    }
}

################################################################################
