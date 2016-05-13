#!/usr/bin/Rscript
library(ggplot2)
library(methods)

# first save pheno.ods into t.in
data=read.table('t.in',header=TRUE)

data$rc_code       = factor(data$rc_code)
data$ascertainment = factor(data$ascertainment)
data$is_female     = factor(data$is_female)

# sample selection

data = subset(data,ldl_x10000>0)

data = subset(data,taking_statins == 0 | is.na(taking_statins))

data = subset(data,stroke_or_tia_diag_age > age_at_sampling |
                   is.na(stroke_or_tia_diag_age) |
                   is.na(age_at_sampling))

#######################################
# age vs ldl
cor(data$age_at_sampling, data$ldl_x10000)
# 0.05979
fit <- lm(ldl_x10000 ~ age_at_sampling, data=data)
summary(fit)
# p-value: 2.621e-15
# very significant

p1 <- ggplot() +
    geom_point(data = data,
               aes(x=age_at_sampling,y=ldl_x10000),
               alpha=0.1) +
    theme_bw() +
    theme(legend.position="none")

#######################################
# pc vs ldl
fit <- lm(ldl_x10000 ~ PC1, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC2, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC3, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC4, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC5, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC6, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC7, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC8, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC9, data=data)
summary(fit)
fit <- lm(ldl_x10000 ~ PC10, data=data)
summary(fit)

# p-value: 0.01788 for PC3,  0.6199 for PC6
# Other PCs are highly correlated.
p2 <- ggplot() +
    geom_point(data = data,
               aes(x=PC2,y=PC1,color=ldl_x10000, fill=ldl_x10000),
               alpha=0.8) +
    ylim(c(0.017,-0.007)) +
    xlim(c(-0.02,0.02)) +
    theme_bw() +
    theme(legend.position="none")

#######################################
# rc vs ldl
bartlett.test(ldl_x10000 ~ rc_code, data=data)
# very significant p-value < 2.2e-16
# heteroscedastic

aov2 <- aov(ldl_x10000 ~ rc_code, data=data)
summary(aov2)
# very significant  p <2e-16

p3 <- ggplot(data, aes(x=rc_code, y=ldl_x10000)) + geom_boxplot() + theme_bw()

#######################################
# gender vs ldl
t.test(ldl_x10000 ~ is_female, data)
# p-value < 2.2e-16
# very significant

p4 <- ggplot(data, aes(x=is_female, y=ldl_x10000)) + geom_boxplot() + theme_bw()

#######################################
# ascertainment vs ldl
aov2 <- aov(ldl_x10000 ~ ascertainment, data=data)
summary(aov2)
# very significant

p5 <- ggplot(data, aes(x=ascertainment, y=ldl_x10000)) +
    geom_boxplot() + theme_bw()

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

png('t1.png',height=11.7, width=8.3,unit='in',res=288)
multiplot(p1, p2,cols=1)
dev.off()


png('t2.png',height=11.7, width=8.3,unit='in',res=288)
multiplot(p3, p4, p5,cols=2)
dev.off()

#######################################
write.table(data,file='t.out',quote=FALSE,row.names=FALSE,sep='\t')


################################################################################
# plot assoc

library(ggplot2)
library(methods)

# head  -n 1 plink.assoc.linear > t.in
# grep ADD plink.assoc.linear | sort -g -k 12 | grep -v NA  >> t.in
# then plot it

data=read.table('t.in',header=TRUE) #  grep ADD plink.assoc.linear

data$minus_log_p = -log10(data$P)

p6 <- ggplot() +
    geom_point(data = data, aes(x=BP,y=minus_log_p,color = R2_1), alpha=0.8) +
    scale_colour_gradient( low="grey", high="red") +
    theme_bw()

p7 <- ggplot() +
    geom_point(data = data, aes(x=BP,y=minus_log_p,color = R2_2), alpha=0.8) +
    scale_colour_gradient( low="grey", high="red") +
    theme_bw()

p8 <- ggplot() +
    geom_point(data = data, aes(x=BP,y=minus_log_p,color = R2_3), alpha=0.8) +
    scale_colour_gradient( low="grey", high="red") +
    theme_bw()

png('t3.png',height=11.7, width=8.3,unit='in',res=288)
multiplot(p6, p7, p8,cols=1)
dev.off()
