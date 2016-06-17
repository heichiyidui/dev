#!/usr/bin/Rscript
library(ggplot2)
library(methods)

args = commandArgs(trailingOnly=TRUE)
ifile_name = args[1]
#  ifile_name = 'rint_p9_0.assoc.linear'
out_man_png_name = paste(ifile_name,'_man.png',sep='')

data=read.table(ifile_name,header=T)
if ('p_1' %in% colnames(data) ) data$P = data$p_1

data=subset(data,!is.na(P))
data$minus_log_p = -log10(data$P)

p1= ggplot(data, aes(x=BP, y=minus_log_p))+
    geom_point(alpha=0.5, size=1.5, shape=19 ) +
    geom_hline( yintercept=-log10(0.0001)) +
    ylab( expression( paste('-log' [10] , italic(P) ) ) ) +
    theme_bw()

png(out_man_png_name,height=6, width=6,unit='in',res=288)
print(p1)
dev.off()

