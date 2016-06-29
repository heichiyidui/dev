#!/usr/bin/Rscript
library(ggplot2)
library(methods)

args = commandArgs(trailingOnly=TRUE)
avm_file_list  = args[1]
# avm_file_list = 'xaa'

avm_files <- scan(avm_file_list, what='')

for (avm_file in avm_files){
    snp_avm  = read.table(avm_file, header = TRUE)
    snp_avm$called = factor(snp_avm$called)

    xtemp <- max(abs(snp_avm$M) )

    snp_cols   = c('0'='#FFCCCC', '1'='#CCCC00', '2'='#00BCFF', '3'='#D4A8FF',
                   '4'='#FF4343', '5'='#949400', '6'='#0058B0', '7'='#A54CFF')
    snp_shapes = c('0' = 16, '1' = 16, '2' = 16, '3' = 16,
                   '4' = 19, '5' = 19, '6' = 19, '7' = 19)
                    # 'AA', 'AB', 'BB' and 'missing'
                    # and highlighted

                    # must be specific, otherwise if one of 'AA' etc is not
                    # presented, 'missing' will be in a different colour.

    p <- ggplot() +
        geom_point(data = subset(snp_avm,called %in% c('0','1','2','3')),
                   aes(x=M, y=A, shape=called, color=called, fill=called),
                   alpha=0.6, size=3) +
        geom_point(data = subset(snp_avm,called %in% c('4','5','6','7')),
                   aes(x=M, y=A, shape=called, color=called, fill=called),
                   alpha=0.9, size=5) +
        scale_colour_manual(values=snp_cols) +
        scale_fill_manual  (values=snp_cols) +
        scale_shape_manual (values=snp_shapes) +
        xlim(c(-xtemp,xtemp)) + xlab(NULL) + ylab(NULL) +
        guides(colour=FALSE,fill=FALSE,shape=FALSE) +
        theme_bw() +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "in"))

    png_file = gsub('.avm','.png',avm_file)
    ggsave(p,file=png_file, dpi=144)
}

