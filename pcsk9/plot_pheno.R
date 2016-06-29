#!/usr/bin/Rscript
library(ggplot2)
library(methods)

# first save pheno.ods into t.in
data=read.table('pheno.csv',header=TRUE)

data$rc      = factor(data$rc)
data$sex     = factor(data$sex)
data$stratum = factor(data$stratum)
data=subset(data,stratum == 6)

#######################################
# age vs ldl
cor(data$age, data$ldl_c)
# 0.0573
fit <- lm(ldl_c ~ age, data=data)
summary(fit)
# p-value: 1.208e-14
# very significant
# not significant for stratum 6? p is 0.03775 there
p1 <- ggplot() +
    geom_point(data = data,
               aes(x=age,y=ldl_c),
               alpha=0.1) +
    theme_bw() +
    theme(legend.position="none")

#######################################
# pc vs ldl

data=read.table('pheno.csv',header=TRUE)
data2=read.table('cov.csv',header=TRUE)
data2$sex <- NULL
data2$FID <- NULL
data2$age <- NULL

data=merge(data,data2,by='IID')
data$rc      = factor(data$rc)
data$sex     = factor(data$sex)
data$stratum = factor(data$stratum)

fit <- lm(ldl_c ~ pc1, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc2, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc3, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc4, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc5, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc6, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc7, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc8, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc9, data=data)
summary(fit)
fit <- lm(ldl_c ~ pc10, data=data)
summary(fit)

# p-value: 0.2403 for pc4, 0.6885 for pc5, 0.1267 for pc10
# Other PCs are highly correlated.

p2 <- ggplot() +
    geom_point(data = data,
               aes(x=pc1,y=pc2,color=ldl_c, fill=ldl_c),
               alpha=0.8) +
    theme_bw() +
    theme(legend.position="none")

#######################################
# rc vs ldl
bartlett.test(ldl_c ~ rc, data=data)
# very significant p-value < 2.2e-16
# heteroscedastic

aov2 <- aov(ldl_c ~ rc, data=data)
summary(aov2)
# very significant  p <2e-16

p3 <- ggplot(data, aes(x=rc, y=ldl_c)) + geom_boxplot() + theme_bw()

#######################################
# gender vs ldl
t.test(ldl_c ~ sex, data)
# p-value = 1.007e-09
# very significant

p4 <- ggplot(data, aes(x=sex, y=ldl_c)) + geom_boxplot() + theme_bw()

#######################################
# ascertainment vs ldl
aov2 <- aov(ldl_c ~ stratum, data=data)
summary(aov2)
# very significant

p5 <- ggplot(data, aes(x=stratum, y=ldl_c)) +
    geom_boxplot() + theme_bw()

################################################################################
# Multiple plot function

# from
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

source('../util/multiplot.R')
################################################################################

png('t1.png',height=11.7, width=8.3,unit='in',res=288)
multiplot(p1, p2,cols=1)
dev.off()


png('t2.png',height=11.7, width=8.3,unit='in',res=288)
multiplot(p3, p4, p5,cols=2)
dev.off()

#######################################
# write.table(data,file='t.out',quote=FALSE,row.names=FALSE,sep='\t')


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
