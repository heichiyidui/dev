#!/usr/bin/Rscript
library(ggplot2)
library(methods)

source('../util/multiplot.R')

data=read.table('t.in',header=T)
data=subset(data,!is.na(meta_p))
data$minus_log_p = -log10(data$meta_p)
data$minus_log_q = -log10(data$meta_q)

p1= ggplot(data, aes(x=BP, y=minus_log_p))+
    geom_point(alpha=0.5, size=1.5, shape=19) +
    geom_hline( yintercept=-log10(0.0001)) +
    ggtitle('Meta-analysis P') +
    theme_bw()

p2= ggplot(data, aes(x=BP, y=minus_log_q))+
    geom_point(alpha=0.5, size=1.5, shape=19) +
    geom_hline( yintercept=-log10(0.05)) +
    ggtitle('Meta-analysis Q') +
    theme_bw()

png('t.png',height=7.5, width=11.7,unit='in',res=288)
multiplot(p1, p2,cols=2)
dev.off()