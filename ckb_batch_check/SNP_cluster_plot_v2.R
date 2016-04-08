#!/usr/bin/Rscript
library(ggplot2)
library(methods)

batch_id  = 'b01'

#######################################
# a funtion to calculate posterior ellipse border coordinates
pos.ellipse <- function(x0, vx, y0, vy, cov)
{
    theta <- 0.5 * atan2(cov * 2, vx - vy)
    sint <- sin(theta)
    cost <- cos(theta)

    a <- 2 * sqrt(vx * cost * cost + vy * sint * sint +
                  cov * 2 * sint * cost)
    b <- 2 * sqrt(vx * sint * sint + vy * cost * cost -
                  cov * 2 * sint * cost)

    np = 100
    alpha <- 2 * pi * (0:np)/np

    sina <- sin(alpha)
    cosa <- cos(alpha)
    x <- x0 + a * cosa * cost - b * sina * sint
    y <- y0 + a * cosa * sint + b * sina * cost

    return (data.frame(x,y))
}

#######################################

ps_file <- paste(batch_id,'.posterior',sep='')

snp_pos <- read.table(ps_file, header = TRUE, as.is =1 )
rownames(snp_pos) <- snp_pos$id
snp_ids <- snp_pos$id

for (snp_id in snp_ids)
{
    el_1 = pos.ellipse(x0  = snp_pos[snp_id,  'x1'],
                       vx  = snp_pos[snp_id, 'vx1'],
                       y0  = snp_pos[snp_id,  'y1'],
                       vy  = snp_pos[snp_id, 'vy1'],
                       cov = snp_pos[snp_id,'cov1'])
    el_2 = pos.ellipse(x0  = snp_pos[snp_id,  'x2'],
                       vx  = snp_pos[snp_id, 'vx2'],
                       y0  = snp_pos[snp_id,  'y2'],
                       vy  = snp_pos[snp_id, 'vy2'],
                       cov = snp_pos[snp_id,'cov2'])
    el_3 = pos.ellipse(x0  = snp_pos[snp_id,  'x3'],
                       vx  = snp_pos[snp_id, 'vx3'],
                       y0  = snp_pos[snp_id,  'y3'],
                       vy  = snp_pos[snp_id, 'vy3'],
                       cov = snp_pos[snp_id,'cov3'])

    avm_file = paste(batch_id,'/',snp_id,'.avm',sep='')
    snp_avm = read.table(avm_file, header = TRUE)
    snp_avm$called = factor(snp_avm$called)

    x3  = snp_pos[snp_id, 'x3']
    vx3 = snp_pos[snp_id,'vx3']
    x1  = snp_pos[snp_id, 'x1']
    vx1 = snp_pos[snp_id,'vx1']

    xtemp <- max(abs(c(x3 + 2*sqrt(vx3), x1 - 2*sqrt(vx1), snp_avm$M) ))

    p <- ggplot() +
        geom_point(data=snp_avm,
                   aes(x=M, y=A, shape=called, color=called, fill=called),
                   alpha=0.5, size=3) +
        scale_colour_manual(values=c('red','#999900','blue','#660066')) +
        scale_fill_manual  (values=c('red','#999900','blue','#660066')) +
        scale_shape_manual (values=c(24,21,25,22) ) +
        geom_path(data=el_1,aes(x=x,y=y),color='blue') +
        geom_path(data=el_2,aes(x=x,y=y),color='#999900') +
        geom_path(data=el_3,aes(x=x,y=y),color='red') +
        xlim(c(-xtemp,xtemp)) + xlab(NULL) + ylab(NULL) +
        guides(colour=FALSE,fill=FALSE,shape=FALSE) +
        theme_bw() +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "in"))

    png_file = paste(batch_id,'/',snp_id,'.png',sep='')
    ggsave(p,file=png_file, dpi=100)
}

