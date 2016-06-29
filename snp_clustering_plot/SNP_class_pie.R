#!/usr/bin/Rscript
library(ggplot2)

snp_cnts = c(521250, 104752, 70838, 1227, 3437, 77710, 2723)
snp_labs = c("PolyHighResolution","NoMinorHom","MonoHighResolution",
             "Hemizygous","OTV","Other","CallRateBelowThreshold")

data <- data.frame(cnts=snp_cnts,type=snp_labs)

col_map = c('PolyHighResolution'    ='#18C45A',
            'NoMinorHom'            ='#5795FF',
            'MonoHighResolution'    ='#18C4B0',
            'Hemizygous'            ='#1874C4',
            'OTV'                   ='#C418C1',
            'Other'                 ='#C41865',
            'CallRateBelowThreshold'='#F36303'
           )

p <- ggplot(data, aes(x = '', y = cnts, fill = type, label=type) ) +
    geom_bar(stat="identity", position = 'fill', width=1 , colour='black')  +
    coord_polar(theta = 'y') +
    scale_fill_manual(values=col_map) +
    xlab('') + ylab('') +
    scale_x_discrete(breaks=NULL) +
    ggtitle("SNP classes") +
    theme_bw()

png('SNP_class_pie.png', height=8.3, width=8.3,unit='in', res=288)
print(p)
dev.off()