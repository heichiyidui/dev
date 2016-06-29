#!/usr/bin/Rscript
library(ggplot2)
library(methods)

data = read.table('manual_chk_res.table',header=TRUE)
# require(gdata)
# data = read.xls('manual_chk_res.xlsx')

data$pass_man_qc[data$pass_man_qc==0] <- 'fail'
data$pass_man_qc[data$pass_man_qc==1] <- 'pass'

snp_cols   = c('fail' = 'red', 'pass' = 'green')

################################################################################
# Multiple plot function
# from
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
source('../util/multiplot.R')

################################################################################
# batch-p vs min-MAF
p1 <- ggplot() +
    geom_point(data = data,
               aes(x=min_maf,y=batch_min_p,color=pass_man_qc,shape=pass_man_qc),
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")

p2 <- ggplot() +
    geom_point(data = data,
               aes(x=min_maf,y=batch_min_p,color=pass_man_qc,shape=pass_man_qc),
               alpha=0.5) +
    xlim(c(0,0.01)) +
    theme_bw() +
    theme(legend.position="none")

png('p1.png', width=8.3, height=11.7, unit='in', res=288)
multiplot(p1, p2,cols=2)
dev.off()

# plate-p vs min-MFA
p1 <- ggplot() +
    geom_point(data = data,
               aes(x=min_maf,y=plate_min_p,color=pass_man_qc,shape=pass_man_qc),
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")

p2 <- ggplot() +
    geom_point(data = data,
               aes(x=min_maf,y=plate_min_p,color=pass_man_qc,shape=pass_man_qc),
               alpha=0.5) +
    xlim(c(0,0.01)) +
    theme_bw() +
    theme(legend.position="none")

png('p2.png', width=8.3, height=11.7, unit='in', res=288)
multiplot(p1, p2,cols=2)
dev.off()


#
# p <- ggplot() +
#     geom_point(data = data,
#                aes(x=min_maf, y=batch_norel_min_p,
#                    color=pass_man_qc, shape=pass_man_qc) ,
#                alpha=0.5) +
#     theme_bw()
# p

# p <- ggplot() +
#     geom_point(data = data,
#                aes(x=batch_min_p,y=plate_min_p,
#                    color=pass_man_qc,shape=pass_man_qc) ,
#                alpha=0.5) +
#     theme_bw()
# p

################################################################################
# batch-p vs max missing call rates
p1 <- ggplot() +
    geom_point(data = data,
               aes(x=max_miss, y=batch_min_p,
                   color=pass_man_qc, shape=pass_man_qc) ,
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")

p2 <- ggplot() +
    geom_point(data = data,
               aes(x=max_miss, y=batch_min_p,
                   color=pass_man_qc, shape=pass_man_qc) ,
               alpha=0.5) +
    xlim(c(0,0.02)) +
    theme_bw() +
    theme(legend.position="none")

png('p3.png', width=8.3, height=11.7, unit='in', res=288)
multiplot(p1, p2,cols=2)
dev.off()

# plate-p vs max missing-call rate

p1 <- ggplot() +
    geom_point(data = data,
               aes(x=max_miss, y=plate_min_p,
                   color=pass_man_qc, shape=pass_man_qc) ,
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")

p2 <- ggplot() +
    geom_point(data = data,
               aes(x=max_miss, y=plate_min_p,
                   color=pass_man_qc, shape=pass_man_qc) ,
               alpha=0.5) +
    xlim(c(0,0.02)) +
    theme_bw() +
    theme(legend.position="none")

png('p4.png', width=8.3, height=11.7, unit='in', res=288)
multiplot(p1, p2,cols=2)
dev.off()

################################################################################
# density

p1 <- ggplot(data, aes(x=batch_min_p, colour=pass_man_qc)) +
    geom_density() + theme_bw() +
    theme(legend.position="none")

p2 <- ggplot(data, aes(x=plate_min_p, colour=pass_man_qc)) +
    geom_density() + theme_bw() +
    theme(legend.position="none")

png('p5.png', width=8.3, height=11.7, unit='in', res=288)
multiplot(p1, p2,cols=1)
dev.off()

p1 <- ggplot(subset(data,batch_min_p!='NA'),
             aes(x=max_miss, colour=pass_man_qc) ) +
    geom_density() +xlim(c(0,0.1)) + theme_bw() +
    theme(legend.position="none")

p2 <- ggplot(subset(data,plate_min_p!='NA'),
             aes(x=max_miss, colour=pass_man_qc) ) +
    geom_density() +xlim(c(0,0.1)) + theme_bw() +
    theme(legend.position="none")
png('p6.png', width=8.3, height=11.7, unit='in', res=288)
multiplot(p1, p2,cols=1)
dev.off()

p1 <- ggplot(subset(data,batch_min_p!='NA'),
             aes(x=min_maf, colour=pass_man_qc) ) +
    geom_density() +xlim(c(0,0.1)) + theme_bw() +
    theme(legend.position="none")

p2 <- ggplot(subset(data,plate_min_p!='NA'),
             aes(x=min_maf, colour=pass_man_qc) ) +
    geom_density() +xlim(c(0,0.1)) + theme_bw() +
    theme(legend.position="none")
png('p7.png', width=8.3, height=11.7, unit='in', res=288)
multiplot(p1, p2,cols=1)
dev.off()

################################################################################
# number of hits

p1 <- ggplot() +
    geom_point(data = data,
               aes(x=batch_num_detected, y=batch_min_p,
                   color=pass_man_qc, shape=pass_man_qc) ,
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")
# SNPs failed two or more batches failed all

p2 <- ggplot() +
    geom_point(data = data,
               aes(x=plate_num_detected, y=plate_min_p,
                   color=pass_man_qc, shape=pass_man_qc) ,
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")
# SNPs failed 4 or more plates failed all

png('p8.png', width=8.3, height=11.7, unit='in', res=288)
multiplot(p1, p2,cols=1)
dev.off()

################################################################################
# heat map

# #ggplot(data,aes(x=max_miss,y=plate_min_p))  + geom_bin2d()

# heatmap_data = read.table('t3.in',header=TRUE)

# png('p9.png', width=8.3, height=11.7, unit='in', res=288)
#     ggplot(heatmap_data,aes(x=max_miss,y=plate_min_p,fill=pass_man_qc)) +
#         geom_tile(colour = "white") +
#         scale_fill_gradient(low = "#CDCDFF", high = "#0000FF") +
#         theme_bw() +        theme(panel.grid = element_blank())

# dev.off()

################################################################################
# plate effect checking

data = read.table('plate_man_qc.table',header=TRUE)

data$pass_qc[data$pass_qc==0] <- 'fail'
data$pass_qc[data$pass_qc==1] <- 'pass'

snp_cols   = c('fail' = 'red', 'pass' = 'green')

p1 <- ggplot() +
    geom_point(data = data,
               aes(x=plate_maf,y=plate_p,color=pass_qc,shape=pass_qc),
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")
# maf doesn't matter

p2 <- ggplot() +
    geom_point(data = data,
               aes(x=plate_maf,y=plate_miss,color=pass_qc,shape=pass_qc),
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")
# missing is the key
p3<- ggplot() +
    geom_point(data = data,
               aes(x=plate_p,y=plate_miss,color=pass_qc,shape=pass_qc),
               alpha=0.5) +
    theme_bw() +
    theme(legend.position="none")
# missing is the key, plate_p somewhat matters


png('p10.png', width=8.3, height=11.7, unit='in', res=288)
    multiplot(p1, p2, p3, cols=1)
dev.off()

pdf('p10.pdf',width=8.3, height = 11.7)
    multiplot(p1, p2, p3, cols=1)
dev.off()

################################################################################
