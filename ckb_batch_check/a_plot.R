#!/usr/bin/Rscript
library(SNPolisher)
library(methods)

snp_file  = 'snp.ls'
batch_id  = 'b01'
call_dir  = 'example_1_buffalo/'

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

to_classify_SNPs = TRUE

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

s = 'ps100'
#for (s in pid)
{
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


SNPolisher:::plot.cluster
function (pid, g, p, col.by = c("called", "reference"),
  type = c("alleleIntensity","AvM"),
  prior = NULL, geno.col, nclus = 3, ...)
{
    col.by <- match.arg(col.by)
    type <- match.arg(type)
    if (col.by == "called") {
        g1 <- g
        g1 <- g1[order(g1[[5]]), ]
        v <- g1$called
    }
    else {
        g1 <- g[g$reference != -9, ]
        v <- g1$reference
    }
    if (col.by == "called") {
        if (is.null(g1$samps.highlight)) {
            g1$bg <- ifelse(v == 0, geno.col[1], ifelse(v ==
                1, geno.col[2], ifelse(v == 2, geno.col[3], ifelse(v ==
                -2, geno.col[5], geno.col[4]))))
        }
        else g1$bg <- ifelse(g1$samps.highlight != "0" & g1$samps.highlight !=
            FALSE, g1$samps.highlight, ifelse(v == 0, geno.col[1],
            ifelse(v == 1, geno.col[2], ifelse(v == 2, geno.col[3],
                ifelse(v == -2, geno.col[5], geno.col[4])))))
    }
    else g1$bg <- ifelse(v == 0, geno.col[1], ifelse(v == 1,
        geno.col[2], ifelse(v == 2, geno.col[3], ifelse(v ==
            -2, geno.col[5], geno.col[4]))))
    g1$col <- ifelse(v == 3, geno.col[7], ifelse(v == 4, geno.col[8],
        "black"))
    g1$pch <- ifelse(v == 0, 24, ifelse(v == 1, 21, ifelse(v ==
        2, 25, ifelse(v == 3, 3, ifelse(v == 4, 4, ifelse(v ==
        -2, 23, 22))))))
    if (type == "alleleIntensity") {
        lim <- c(0, max(g$a_signal, g$b_signal))
        plot(0, 0, type = "n", xlim = lim, ylim = lim, xlab = "A",
            ylab = "B", main = paste(pid, "\nalleleSignal -",
                col.by, "genotypes", sep = " "), ...)
        points(g1$a_signal, g1$b_signal, col = g1$col, pch = g1$pch,
            bg = g1$bg)
    }
    else {
        a <- log2(g1$a_signal)
        b <- log2(g1$b_signal)
        A <- (a + b)/2
        M <- a - b
        xtemp <- max(abs(c(p$x[3] + 2 * sqrt(p$vx[3]), p$x[1] -
            2 * sqrt(p$vx[1]), M)))
        ylim <- range(c(p$y + 2 * sqrt(p$vy), p$y - 2 * sqrt(p$vy),
            A))
        plot(M, A, xlim = c(-xtemp, xtemp), ylim = ylim, xlab = "A-B",
            ylab = "(A+B)/2", main = paste(pid, "\nAvM -", col.by,
                "genotypes", sep = " "), col = g1$col, pch = g1$pch,
            bg = g1$bg, cex.main = 1.2, ...)
        if (col.by == "called")
            sapply(1:nclus, function(j) {
                plot.cov(p$vx[j], p$vy[j], p$cov[j], p$x[j],
                  p$y[j], lwd = 2, col = "#AAAA44")
            })
        if (!is.null(prior) & col.by == "called")
            sapply(1:nclus, function(j) {
                with(prior, plot.cov(vx[j], vy[j], cov[j], x[j],
                  y[j], lt = 2, lwd = 2, col = "#44AAAA"))
            })
    }
}


s
[1] "ps100"
g
    called  a_signal  b_signal    sample samps.highlight reference
1        1  3203.898  3883.094   Sample1           FALSE        -1
2        1  4315.499  5930.525   Sample2           FALSE        -1
3        1  3879.043  4475.433   Sample3           FALSE        -1
4        0 10014.346  1720.226   Sample4           FALSE        -1
5        0 10355.236  1545.243   Sample5           FALSE        -1
6        1  3582.515  4891.744   Sample6           FALSE        -1
7        2  1326.958 12194.912   Sample7           FALSE        -1
8        1  4314.096  6990.968   Sample8           FALSE        -1
9        1  3136.829  4636.816   Sample9           FALSE        -1
10       0  6564.856  1848.515  Sample10           FALSE        -1
11       2  1430.507  8683.613  Sample11           FALSE        -1
12       1  3340.885  4289.189  Sample12           FALSE        -1
13       1  4568.688  5198.452  Sample13           FALSE        -1
14       0  8920.109  1751.875  Sample14           FALSE        -1
15       1  4174.439  5019.193  Sample15           FALSE        -1
16       2  1673.044 10683.125  Sample16           FALSE        -1
17       0  8174.149  1624.647  Sample17           FALSE        -1
18       0  8928.374  1743.033  Sample18           FALSE        -1
19       1  4402.606  6025.956  Sample19           FALSE        -1
20       2  1716.718 16286.880  Sample20           FALSE        -1
21       1  3543.287  4392.708  Sample21           FALSE        -1
22       0  7591.509  1644.983  Sample22           FALSE        -1
23       1  3753.720  4506.610  Sample23           FALSE        -1
24       2  1194.887  5421.055  Sample24           FALSE        -1
25       0  9548.345  1668.547  Sample25           FALSE        -1
26       2  1288.948 10812.393  Sample26           FALSE        -1
27      -1  1938.049  6691.132  Sample27           FALSE        -1

p
          x        vx       y        vy         cov
1 -2.884160 0.0595628 11.8976 0.0375299 -0.03130760
2 -0.399992 0.0595628 12.1704 0.0375299  0.00626921
3  2.234520 0.0595628 11.7516 0.0375299  0.01684390



SNPolisher:::plot.cov
function (vx, vy, cov, x = 0, y = 0, ...)
{
    theta <- 0.5 * atan2(cov * 2, vx - vy)
    sint <- sin(theta)
    cost <- cos(theta)
    yz.ellipse(2 * sqrt(vx * cost * cost + vy * sint * sint +
        cov * 2 * sint * cost), 2 * sqrt(vx * sint * sint + vy *
        cost * cost - cov * 2 * sint * cost), theta, x, y, ...)
}

SNPolisher:::yz.ellipse
function (a = 1, b = 1, theta = 0, x0 = 0, y0 = 0, np = 100,
    ...)
{
    alpha <- 2 * pi * (0:np)/np
    sint <- sin(theta)
    cost <- cos(theta)
    sina <- sin(alpha)
    cosa <- cos(alpha)
    x <- x0 + a * cosa * cost - b * sina * sint
    y <- y0 + a * cosa * sint + b * sina * cost
    lines(x, y, type = "l", ...)
}


library(ggplot2)

dat = data.frame(x,y)
ggplot() +
    geom_path(aes(x,y))