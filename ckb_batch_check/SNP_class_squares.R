#!/usr/bin/Rscript
library(ggplot2)

options(stringsAsFactors = FALSE)

#######################################
# read classifications
snps = read.table('snp.ls',header=FALSE)
snps = c(snps[,1])

snp_classes <- data.frame(snps)
colnames(snp_classes) <- c('probeset_id')

batches = c('b01','b02','b03','b04','b05','b06','b07')

for (batch in batches)
{
    data=read.table(paste(batch,'/Ps.performance.txt',sep=''),header=TRUE)
    sub_data = data[data$probeset_id %in% snps,]
    sub_data = subset(sub_data,select = c(probeset_id,ConversionType))
    snp_classes = merge(snp_classes,sub_data,by='probeset_id')

    colnames(snp_classes)[
        which(names(snp_classes) == "ConversionType")] <-
            paste(batch,'class',sep='_')
}

#######################################
# plot the SNPs

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

# snp = 'AX-100004419'
for (snp in snps)
{
    s_cls = snp_classes[snp_classes$probeset_id==snp,]
    s_cls$probeset_id <- NULL
    s_cls = c(t(s_cls))    # SNP classes
    fills = col_map[s_cls] # filling colours

    png(paste('class_png/',snp,'_class.png',sep=''))

    p <- ggplot(data=rec_location, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2) ) +
        geom_rect(fill=fills) +
        xlab('') + ylab('') +
        scale_x_continuous(breaks=NULL) +
        scale_y_continuous(breaks=NULL) +
        theme_bw()

    print(p)
    dev.off()
}

# the end
################################################################################

