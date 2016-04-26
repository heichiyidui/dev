#!/usr/bin/Rscript
library(ggplot2)
library(methods)

data = read.table('t.in',header=TRUE)
# require(gdata)
# data = read.xls('manual_chk_res.xlsx')

data$pass_man_qc[data$pass_man_qc==0] <- 'fail'
data$pass_man_qc[data$pass_man_qc==1] <- 'pass'

snp_cols   = c('fail' = 'red', 'pass' = 'green')

################################################################################
# Multiple plot function
# from
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# ggplot objects can be passed in ..., or to plotlist
# (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
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

png('p1.png',width=595,height=842)
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

png('p2.png',width=595,height=842)
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

png('p3.png',width=595,height=842)
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

png('p4.png',width=595,height=842)
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

png('p5.png',width=595,height=842)
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
png('p6.png',width=595,height=842)
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
png('p7.png',width=595,height=842)
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

png('p8.png',width=595,height=842)
multiplot(p1, p2,cols=1)
dev.off()

################################################################################
# heat map

#ggplot(data,aes(x=max_miss,y=plate_min_p))  + geom_bin2d()

heatmap_data = read.table('t3.in',header=TRUE)

png('p9.png',width=595,height=842)
    ggplot(heatmap_data,aes(x=max_miss,y=plate_min_p,fill=pass_man_qc)) +
        geom_tile(colour = "white") +
        scale_fill_gradient(low = "#CDCDFF", high = "#0000FF") +
        theme_bw() +        theme(panel.grid = element_blank())

dev.off()
################################################################################