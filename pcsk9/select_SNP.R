#!/usr/bin/Rscript
library(ggplot2)
library(methods)

# data1=read.table('p9.dosage',header=T)
# data2=read.table('pheno.csv',header=T)

# data = merge(data1,data2,by='IID')

# data_features = subset(data,select=c(AX.83389438,AX.39912161,AX.31642001,
#                                      AX.31642169,AX.11541856,AX.11576926))
# data_target = t(subset(data,select=c(rint_ldl_c)))

# library(caret)

# data_features<- scale(data_features, center=TRUE, scale=TRUE);

# inTrain <- createDataPartition(data_target, p = 3/4, list = FALSE);
# # 16702 in training

# #Create the Training Dataset for Descriptors
# trainDescr <- data_features[inTrain,];

# # Create the Testing dataset for Descriptors
# testDescr <- data_features[-inTrain,];

# trainTarget <- data_target[inTrain];
# testTarget  <- data_target[-inTrain];


# descrCorr <- cor(trainDescr,use = "na.or.complete" );
# highCorr <- findCorrelation(descrCorr, 0.70);
# trainDescr <- trainDescr[, -highCorr];
# testDescr <- testDescr[, -highCorr];
# # Here, we can included a correlation matrix analysis to remove the redundant
# # features before the backwards selection
# # AX.31642169 removed!

# svmProfile <- rfe(x=trainDescr, y = trainTarget, sizes = c(1:5),
#                   rfeControl= rfeControl(functions = caretFuncs,number = 2),
#                   method = "svmRadial",fit = FALSE);
# # caret function: the rfe is the backwards selection,
# # c is the possible sizes of the features sets,
# # and method the optimization method is a support vector machine.

# What the heck? It gave a single SNP AX-83389438 after hours.
################################################################################
# The SVM training is running forever
# stepAIC?

data1=read.table('p9.dosage',header=T)
data2=read.table('pheno.csv',header=T)

data2$rc      = as.factor(data2$rc)
data2$stratum = as.factor(data2$stratum)

data = merge(data1,data2,by='IID')
data = data[rowSums(is.na(data)) == 0,]
data = subset(data,stratum!=6)

library(MASS)
fit <-lm(ldl_c ~ rc + sex + age + stratum +
                 AX.83389438 + AX.39912161 + AX.31642001 +
                 AX.31642169 + AX.11541856 + AX.11576926,
         data=data)

# fit <-lm(ldl_c ~ rc + sex + age + stratum +
#                  AX.83389438+AX.39912161+AX.31642001+AX.11576926+AX.31642169+
#                  AX.31641677+AX.11541856+AX.31641243+AX.39912159+AX.51209582+
#                  AX.31642257+AX.11629248+AX.11307124+AX.12935655+AX.11150762+
#                  AX.31648117+AX.12932322+AX.83187812+AX.50958112+AX.11447635,
#          data=data)

step <- stepAIC(fit, direction = "backward")
# or step <- stepAIC(fit, direction = "both")
step$anova


data = merge(data1,data2,by='IID')
data = data[rowSums(is.na(data)) == 0,]

fit <- lm(rint_ldl_c ~ AX.83389438 + AX.39912161 + AX.31642001 +
                       AX.31642169 + AX.11541856 + AX.11576926,
          data=data)

step <- stepAIC(fit, direction = "backward")
step$anova
