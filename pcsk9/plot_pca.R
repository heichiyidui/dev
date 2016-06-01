#!/usr/bin/Rscript
library(ggplot2)
library(methods)

data=read.table('t.in',header=T)

data$rc = as.factor(data$rc)

rc_cols = c('12'='green2', '16'='goldenrod1', '26'='magenta', '36'='cyan2',
            '46'='pink1',  '52'='red',        '58'='gold4',   '68'='dimgrey',
            '78'='blue',   '88'='red4')

rc_names = c('12'='Qingdao',  '16'='Harbin',  '26'='Haikou', '36'='Suzhou',
             '46'='Liuzhou',  '52'='Sichuan', '58'='Gansu',  '68'='Henan',
             '78'='Zhejiang', '88'='Hunan')

source('../util/multiplot.R')

# family members against the others
p1 = ggplot(data, aes(x=pc1, y=pc2, color=is_fam)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    theme_bw() +
    theme(legend.position="none")

p2 = ggplot(data, aes(x=pc3, y=pc4, color=is_fam)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    theme_bw() +
    theme(legend.position="none")

p3 = ggplot(data, aes(x=pc5, y=pc6, color=is_fam)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    theme_bw() +
    theme(legend.position="none")

p4 = ggplot(data, aes(x=pc7, y=pc8, color=is_fam)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    theme_bw() +
    theme(legend.position="none")

p5 = ggplot(data, aes(x=pc9, y=pc10, color=is_fam)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    theme_bw() +
    theme(legend.position="none")

png('pca_fam.png',height=11.7, width=7.8,unit='in',res=288)
multiplot(p1, p2, p3, p4, p5, cols=2)
dev.off()


# rc vs pcs
p1 = ggplot(data, aes(x=pc1, y=pc2, color=rc)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    scale_colour_manual(values=rc_cols) +
    theme_bw() +
    theme(legend.position="none")

p2 = ggplot(data, aes(x=pc3, y=pc4, color=rc)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    scale_colour_manual(values=rc_cols) +
    theme_bw() +
    theme(legend.position="none")

p3 = ggplot(data, aes(x=pc5, y=pc6, color=rc)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    scale_colour_manual(values=rc_cols) +
    theme_bw() +
    theme(legend.position="none")

p4 = ggplot(data, aes(x=pc7, y=pc8, color=rc)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    scale_colour_manual(values=rc_cols) +
    theme_bw() +
    theme(legend.position="none")

p5 = ggplot(data, aes(x=pc9, y=pc10, color=rc)) +
    geom_point(alpha=0.5, size=1.5, shape=19) +
    scale_colour_manual(values=rc_cols) +
    theme_bw() +
    theme(legend.position="none")


png('pca_rc.png',height=11.7, width=7.8,unit='in',res=288)
multiplot(p1, p2, p3, p4, p5, cols=2)
dev.off()

