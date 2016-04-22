#!/usr/bin/Rscript
library(ggplot2)
library(methods)

data = read.table('t.in',header=TRUE)

data$pass_man_qc[data$pass_man_qc==0] <- 'fail'
data$pass_man_qc[data$pass_man_qc==1] <- 'pass'

snp_cols   = c('fail' = 'red', 'pass' = 'green')

p <- ggplot() +
    geom_point(data = data,
               aes(x=min_maf, y=batch_min_p, color=pass_man_qc) ,
               alpha=0.5) +
    theme_bw()
p


p <- ggplot() +
    geom_point(data = data,
               aes(x=min_maf, y=plate_min_p, color=pass_man_qc) ,
               alpha=0.5) +
    theme_bw()
p

# AX-105017241 batch_norel_min_p 8.35 -> NA

p <- ggplot() +
    geom_point(data = data,
               aes(x=min_maf, y=batch_norel_min_p, color=pass_man_qc) ,
               alpha=0.5) +
    theme_bw()
p


p <- ggplot() +
    geom_point(data = data,
               aes(x=batch_min_p, y=plate_min_p, color=pass_man_qc) ,
               alpha=0.5) +
    theme_bw()
p


p <- ggplot() +
    geom_point(data = data,
               aes(x=max_miss, y=plate_min_p, color=pass_man_qc) ,
               alpha=0.5) +
    theme_bw()
p

p <- ggplot() +
    geom_point(data = data,
               aes(x=max_miss, y=batch_min_p, color=pass_man_qc) ,
               alpha=0.5) +
    theme_bw()
p

