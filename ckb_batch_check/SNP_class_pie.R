#!/usr/bin/Rscript
library(ggplot2)


col_map = vector(mode='list',length=7)
names(col_map) = c("PolyHighResolution", "NoMinorHom", "MonoHighResolution",
                   "Hemizygous", "OTV", "Other", "CallRateBelowThreshold")

col_map[['PolyHighResolution']]     <- '#18C45A'
col_map[['NoMinorHom']]             <- '#5795FF'
col_map[['MonoHighResolution']]     <- '#18C4B0'
col_map[['Hemizygous']]             <- '#1874C4'
col_map[['OTV']]                    <- '#C418C1'
col_map[['Other']]                  <- '#C41865'
col_map[['CallRateBelowThreshold']] <- '#F36303'


snp_cnts = c(521250, 104752, 70838, 1227, 3437, 77710, 2723)
snp_labs = c("PolyHighResolution","NoMinorHom","MonoHighResolution",
             "Hemizygous","OTV","Other","CallRateBelowThreshold")

dat <- data.frame(cnts=snp_cnts,type=snp_labs)

cols = c('#F36303','#5795FF','#18C4B0','#1874C4','#C418C1','#C41865','#18C45A')

png('SNP_class_pie.png')

p <- ggplot(dat, aes(x = '', y = cnts, fill = type, label=type ) ) +
    geom_bar(stat="identity", position = 'fill', width=1 , colour='black')  +
    coord_polar(theta = 'y') +
    scale_fill_manual(values=cols) +
    xlab('') + ylab('') +
    scale_x_discrete(breaks=NULL) +
    ggtitle("SNP classes") +
    theme_bw()
print(p)

dev.off()