#!/usr/bin/Rscript
library(methods)

# #######################################
# # to cluster AA_AA contacts using pt120 and plot them
# library('fastcluster')

# cont_labels = c(
# 'AA','AR','AN','AD','AC','AQ','AE','AG','AH','AI','AL','AK','AM','AF','AP','AS',
# 'AT','AW','AY','AV',
# 'RA','RR','RN','RD','RC','RQ','RE','RG','RH','RI','RL','RK','RM','RF','RP','RS',
# 'RT','RW','RY','RV',
# 'NA','NR','NN','ND','NC','NQ','NE','NG','NH','NI','NL','NK','NM','NF','NP','NS',
# 'NT','NW','NY','NV',
# 'DA','DR','DN','DD','DC','DQ','DE','DG','DH','DI','DL','DK','DM','DF','DP','DS',
# 'DT','DW','DY','DV',
# 'CA','CR','CN','CD','CC','CQ','CE','CG','CH','CI','CL','CK','CM','CF','CP','CS',
# 'CT','CW','CY','CV',
# 'QA','QR','QN','QD','QC','QQ','QE','QG','QH','QI','QL','QK','QM','QF','QP','QS',
# 'QT','QW','QY','QV',
# 'EA','ER','EN','ED','EC','EQ','EE','EG','EH','EI','EL','EK','EM','EF','EP','ES',
# 'ET','EW','EY','EV',
# 'GA','GR','GN','GD','GC','GQ','GE','GG','GH','GI','GL','GK','GM','GF','GP','GS',
# 'GT','GW','GY','GV',
# 'HA','HR','HN','HD','HC','HQ','HE','HG','HH','HI','HL','HK','HM','HF','HP','HS',
# 'HT','HW','HY','HV',
# 'IA','IR','IN','ID','IC','IQ','IE','IG','IH','II','IL','IK','IM','IF','IP','IS',
# 'IT','IW','IY','IV',
# 'LA','LR','LN','LD','LC','LQ','LE','LG','LH','LI','LL','LK','LM','LF','LP','LS',
# 'LT','LW','LY','LV',
# 'KA','KR','KN','KD','KC','KQ','KE','KG','KH','KI','KL','KK','KM','KF','KP','KS',
# 'KT','KW','KY','KV',
# 'MA','MR','MN','MD','MC','MQ','ME','MG','MH','MI','ML','MK','MM','MF','MP','MS',
# 'MT','MW','MY','MV',
# 'FA','FR','FN','FD','FC','FQ','FE','FG','FH','FI','FL','FK','FM','FF','FP','FS',
# 'FT','FW','FY','FV',
# 'PA','PR','PN','PD','PC','PQ','PE','PG','PH','PI','PL','PK','PM','PF','PP','PS',
# 'PT','PW','PY','PV',
# 'SA','SR','SN','SD','SC','SQ','SE','SG','SH','SI','SL','SK','SM','SF','SP','SS',
# 'ST','SW','SY','SV',
# 'TA','TR','TN','TD','TC','TQ','TE','TG','TH','TI','TL','TK','TM','TF','TP','TS',
# 'TT','TW','TY','TV',
# 'WA','WR','WN','WD','WC','WQ','WE','WG','WH','WI','WL','WK','WM','WF','WP','WS',
# 'WT','WW','WY','WV',
# 'YA','YR','YN','YD','YC','YQ','YE','YG','YH','YI','YL','YK','YM','YF','YP','YS',
# 'YT','YW','YY','YV',
# 'VA','VR','VN','VD','VC','VQ','VE','VG','VH','VI','VL','VK','VM','VF','VP','VS',
# 'VT','VW','VY','VV')

# pt120=as.matrix(read.table('050',header=F))
# # it's actually the sum of sum_mat 040~060

# rownames(pt120) = cont_labels
# colnames(pt120) = cont_labels

# dis=as.dist(-log(pt120),diag=F,upper=F)

# tree = hclust(dis, method='ward.D2' )

# pdf('cont_groups_1.pdf', width = 28)
# par(cex=0.5)
# plot(tree)

# clu = rect.hclust(tree,k=20)
# dev.off()

#######################################
# Contacts are clustered into 20 groups, according to their exchangeablity
# in the sum_mutation matrix.

#######################################
# using different colours

library(ggplot2)

data=read.table('t.in',header=TRUE)
data$group = as.factor(data$group)

aa_labels = c('A','R','N','D','C','Q','E','G','H','I',
              'L','K','M','F','P','S','T','W','Y','V')

# roughly 4 groups
g1=subset(data,group %in% c( 0, 1, 2, 3))
g2=subset(data,group %in% c( 4, 5, 6, 7, 8, 9,10))
g3=subset(data,group %in% c(11,12,13,14,15))
g4=subset(data,group %in% c(16,17,18,19))

p1 <- ggplot(g1, aes(x=r2, y=r1, fill=group, width=0.9, height=0.9)) +
    geom_tile() +
    xlab('second residue') +
    ylab('first residue') +
    scale_x_discrete(limits=aa_labels,drop=F) +
    scale_y_discrete(limits=aa_labels,drop=F) +
    ggtitle("Group 1") +
    theme_bw() +
    theme(legend.position="none")

p2 <- ggplot(g2, aes(x=r2, y=r1, fill=group, width=0.9, height=0.9 )) +
    geom_tile() +
    xlab('second residue') +
    ylab('first residue') +
    scale_x_discrete(limits=aa_labels,drop=F) +
    scale_y_discrete(limits=aa_labels,drop=F) +
    ggtitle("Group 2") +
    theme_bw() +
    theme(legend.position="none")

p3 <- ggplot(g3, aes(x=r2, y=r1, fill=group, width=0.9, height=0.9 )) +
    geom_tile() +
    xlab('second residue') +
    ylab('first residue') +
    scale_x_discrete(limits=aa_labels,drop=F) +
    scale_y_discrete(limits=aa_labels,drop=F) +
    ggtitle("Group 3") +
    theme_bw() +
    theme(legend.position="none")

p4 <- ggplot(g4, aes(x=r2, y=r1, fill=group, width=0.9, height=0.9 )) +
    geom_tile() +
    xlab('second residue') +
    ylab('first residue') +
    scale_x_discrete(limits=aa_labels,drop=F) +
    scale_y_discrete(limits=aa_labels,drop=F) +
    ggtitle("Group 4") +
    theme_bw() +
    theme(legend.position="none")

################################################################################
# Multiple plot function
# from
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
source('../util/multiplot.R')
################################################################################

pdf('cont_groups_2.pdf',height=8.3, width=8.3)
multiplot(p1, p3, p2,p4,cols=2)
dev.off()
