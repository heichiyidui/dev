#!/usr/bin/Rscript
library(methods)
library('fastcluster')

#######################################
# to cluster AA_AA contacts and plot them
# from pt120

cont_labels = c(
'AA','AR','AN','AD','AC','AQ','AE','AG','AH','AI','AL','AK','AM','AF','AP','AS','AT','AW','AY','AV',
'RA','RR','RN','RD','RC','RQ','RE','RG','RH','RI','RL','RK','RM','RF','RP','RS','RT','RW','RY','RV',
'NA','NR','NN','ND','NC','NQ','NE','NG','NH','NI','NL','NK','NM','NF','NP','NS','NT','NW','NY','NV',
'DA','DR','DN','DD','DC','DQ','DE','DG','DH','DI','DL','DK','DM','DF','DP','DS','DT','DW','DY','DV',
'CA','CR','CN','CD','CC','CQ','CE','CG','CH','CI','CL','CK','CM','CF','CP','CS','CT','CW','CY','CV',
'QA','QR','QN','QD','QC','QQ','QE','QG','QH','QI','QL','QK','QM','QF','QP','QS','QT','QW','QY','QV',
'EA','ER','EN','ED','EC','EQ','EE','EG','EH','EI','EL','EK','EM','EF','EP','ES','ET','EW','EY','EV',
'GA','GR','GN','GD','GC','GQ','GE','GG','GH','GI','GL','GK','GM','GF','GP','GS','GT','GW','GY','GV',
'HA','HR','HN','HD','HC','HQ','HE','HG','HH','HI','HL','HK','HM','HF','HP','HS','HT','HW','HY','HV',
'IA','IR','IN','ID','IC','IQ','IE','IG','IH','II','IL','IK','IM','IF','IP','IS','IT','IW','IY','IV',
'LA','LR','LN','LD','LC','LQ','LE','LG','LH','LI','LL','LK','LM','LF','LP','LS','LT','LW','LY','LV',
'KA','KR','KN','KD','KC','KQ','KE','KG','KH','KI','KL','KK','KM','KF','KP','KS','KT','KW','KY','KV',
'MA','MR','MN','MD','MC','MQ','ME','MG','MH','MI','ML','MK','MM','MF','MP','MS','MT','MW','MY','MV',
'FA','FR','FN','FD','FC','FQ','FE','FG','FH','FI','FL','FK','FM','FF','FP','FS','FT','FW','FY','FV',
'PA','PR','PN','PD','PC','PQ','PE','PG','PH','PI','PL','PK','PM','PF','PP','PS','PT','PW','PY','PV',
'SA','SR','SN','SD','SC','SQ','SE','SG','SH','SI','SL','SK','SM','SF','SP','SS','ST','SW','SY','SV',
'TA','TR','TN','TD','TC','TQ','TE','TG','TH','TI','TL','TK','TM','TF','TP','TS','TT','TW','TY','TV',
'WA','WR','WN','WD','WC','WQ','WE','WG','WH','WI','WL','WK','WM','WF','WP','WS','WT','WW','WY','WV',
'YA','YR','YN','YD','YC','YQ','YE','YG','YH','YI','YL','YK','YM','YF','YP','YS','YT','YW','YY','YV',
'VA','VR','VN','VD','VC','VQ','VE','VG','VH','VI','VL','VK','VM','VF','VP','VS','VT','VW','VY','VV')

pt120=as.matrix(read.table('050',header=F))
# it's actually the sum of sum_mat 040~060

rownames(pt120) = cont_labels
colnames(pt120) = cont_labels

dis=as.dist(-log(pt120),diag=F,upper=F)

tree = hclust(dis, method='ward.D2' )

pdf('t.pdf', width = 28)
par(cex=0.5)
plot(tree)

clu = rect.hclust(tree,k=20)
dev.off()

#######################################
# plot the grouping now

aa_labels = c('A','R','N','D','C','Q','E','G','H','I',
              'L','K','M','F','P','S','T','W','Y','V')
data=as.matrix(read.table('t.in',header=F))

rownames(data) = aa_labels
colnames(data) = aa_labels

heatmap(data)


#######################################
# using different colours, change row and column orders

library(ggplot2)

data=read.table('t.in',header=TRUE)
data$group = as.factor(data$group)

row_labels = c('W','F','Y','P','V','M','I','L','C','G',
               'T','A','S','H','D','N','E','Q','R','K')
col_labels = c('C','T','A','S','V','M','I','L','W','F',
               'Y','P','H','N','Q','E','D','G','R','K')

group_colors = c("#F7FCF5","#F4FAF1","#F1F9EE","#EEF8EA","#C8E9C1",
				 "#C2E7BB","#BCE4B5","#B5E1AF","#AFDFA9","#A9DCA3",
				 "#A3D99D","#4EB163","#46AD5F","#3FA95B","#3AA357",
				 "#359E53","#006528","#005E26","#005823","#005120")
# roughly 4 groups

ggplot(data, aes(x, y, fill = group )) +
    xlab('the second residue') +
    ylab('the first residue') +
    scale_x_discrete(limits=col_labels) +
    scale_y_discrete(limits=row_labels) +
    geom_tile() +
    theme_bw() +
    scale_fill_manual(values = group_colors) +
    theme(legend.position="none")


