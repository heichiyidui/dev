################################################################################
# 1. Ancestries PCA plot 

#######################################
# read table

ance_pca_data= read.table('ancestries.in',as.is=3,header=TRUE)

#######################################
# simple plot 
plot(ance_pca_data$PC1, ance_pca_data$PC2, pch=19, 
     col = factor(ance_pca_data$cohort) ,
     xlab = 'PC1' , ylab = 'PC2', main = 'Ancestries PCA' )
     
# add legend, which is of wrong locations anyway
legend(0,0,unique(ance_pca_data$cohort),
       col=1:length(ance_pca_data$cohort),pch=19)

#######################################
# ggplot2 plotting

library(ggplot2)

ance_pca_data= read.table('ancestries.in',as.is=3,header=TRUE)

# simple plot

qplot(PC1, PC2, data=ance_pca_data, color=cohort,
      main = "Ancestries PCA", 
      size = I(4), alpha=I(0.5) ) + theme_bw()
# too many cohorts

# Changed cohort ids in the table 'ancestries.in' to continent ids and 'CKB'.
# put 'EAS' to the end of the file, behind 'CKB'

ance_pca_data= read.table('ancestries_cont.in',as.is=4,header=TRUE)

# head of file ancestries_cont.in
# PC1 PC2 PC3 cohort
# -0.0200812 -0.026702 0.00952348 EUR
# -0.0205353 -0.0268003 0.00911413 EUR
# -0.0204138 -0.027077 0.0096619 EUR
# ...

pdf('ancestries_pca.pdf',width=8.5)
    qplot(PC1, PC2, data=ance_pca_data, color=cohort,
          ylim = c(0.03,-0.03) ,
          main = "Ancestries PCA", 
          size = I(2), alpha=I(0.4), shape= I(16) ) + 
    scale_colour_manual(
        name = "",
        values=c('red','yellow4','black','cyan2','blue','magenta')) + 
    theme_bw()

dev.off()

pdf('ancestries_pca_pc2vs3.pdf',width=8.5)
    qplot(PC3, PC2, data=ance_pca_data, color=cohort,
          main = "Ancestries PCA", 
          size = I(2), alpha=I(0.4), shape= I(16) ) + 
    scale_colour_manual(
        name = "",
        values=c('red','yellow4','black','cyan2','blue','magenta')) + 
    theme_bw()

dev.off()

################################################################################
# 2. CKB region PCA plots

library(ggplot2)

# Randomly shuffled and selected 7000 subjects for the plots. 
# The orginal file was saved as 'ckb_pca.in.bak' 
data = read.table('ckb_pca.in',header=TRUE,as.is=1)

# simple plot 
pdf('ckb_pca.pdf',width=8.5)
    qplot(PC2, PC1, data = data, color=cohort,
          ylim = c(0.02,-0.01) ,
          main = "CKB PCA", 
          size = I(1.5), alpha=I(0.4), shape= I(16) ) +
    scale_colour_manual(
        name = "",
        values=c('gold4','magenta','goldenrod1','dimgrey','red4',
                 'pink1','green2','red','cyan2','blue')) + 
    theme_bw()

dev.off()

###################
# can be done with ggplot too
# p = ggplot(data, aes(x=PC2, y=PC1, color=cohort))  + 
#     geom_point(alpha=0.4) +
#     ylim( c(0.02,-0.01) ) + theme_bw() 
# p

#######################################
# ggplot with manually chosen colours

# rural cohorts
dat_rural = data[data$cohort %in% 
                 c('Gansu','Henan','Hunan','Sichuan','Zhejiang'),]

pdf('ckb_rural_pca.pdf',width=8.5)
    qplot(PC2, PC1, data = dat_rural, color=cohort,
          ylim = c(0.01,-0.01) ,
          main = "CKB_rural PCA", 
          size = I(1.5), alpha=I(0.4), shape= I(16) ) + 
    scale_colour_manual(name = "",
                        values = c("gold4","dimgrey", "red4","red","blue")) +
    theme_bw()

dev.off()

# urban cohorts
dat_urban = data[data$cohort %in%
                 c('Haikou','Harbin','Liuzhou','Qingdao','Suzhou'),]

pdf('ckb_urban_pca.pdf',width=8.5)
    qplot(PC2, PC1, data = dat_urban, color=cohort,
          ylim = c(0.02,-0.01) ,
          main = "CKB_urban PCA", 
          size = I(1.5), alpha=I(0.4), shape= I(16) ) + 
    scale_colour_manual(name = "",
                        values = c("magenta","goldenrod1","pink1",
                                   "green2","cyan2")) +
    # to remove the grid
    # theme(panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank()) +
    theme_bw() 

dev.off()

################################################################################
# 3. plotting PCA outliers?

qplot(PC8,PC9,data=data) + theme_bw()

################################################################################