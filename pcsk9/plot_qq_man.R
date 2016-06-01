#!/usr/bin/Rscript
library(ggplot2)
library(methods)
source('../util/qqman.R')

args = commandArgs(trailingOnly=TRUE)
ifile_name = args[1]

data=read.table(ifile_name,header=T)
data=subset(data,!is.na(P))

# data=subset(data,P > 1.0e-08)

cbPalette  <- c("#666666", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

out_man_png_name = paste(ifile_name,'_man.png',sep='')
png(out_man_png_name,height=8.3, width=11.7, unit='in', res=288)
manhattan(data,col=cbPalette)
dev.off()

out_qq_png_name = paste(ifile_name,'_qq.png',seq='')
png(out_qq_png_name,res=288)
qq(data$P)
dev.off()


