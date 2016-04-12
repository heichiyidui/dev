#!/usr/bin/Rscript
library(ggplot2)

#######################################
# read class data

data=read.table('t.in', header=FALSE)

colnames(data)=c('snp_id','b01','b02','b03','b04','b05','b06','b07')
rownames(data)=data$snp_id

#######################################
# class colour map

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

rec_location = data.frame(
        x1 = c(0.36,0.69,0.03,0.36,0.69,0.03,0.36),
        x2 = c(0.63,0.96,0.30,0.63,0.96,0.30,0.63),
        y1 = c(0.69,0.69,0.36,0.36,0.36,0.03,0.03),
        y2 = c(0.96,0.96,0.63,0.63,0.63,0.30,0.30)
    )

#######################################
# plot the SNPs

for ( snp_id in data$snp_id )
{
    s_cls = data[snp_id,]
    s_cls$snp_id <- NULL
    s_cls = c(t(s_cls))    # SNP classes
    fills = col_map[s_cls] # filling colours

    png(paste('class_png/',snp_id,'.png',sep=''))

    p <- ggplot(data=rec_location, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2) ) +
        geom_rect(fill=fills) +
        xlab('') + ylab('') +
        scale_x_continuous(breaks=NULL) +
        scale_y_continuous(breaks=NULL) +
        theme_bw() +
        theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "in"))

    png_file = paste('class_png/',snp_id,'.png',sep='')
    ggsave(p,file=png_file, dpi=100)
}

