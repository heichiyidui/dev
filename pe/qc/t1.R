IMISS=read.table("pe9_missing.imiss", header=T, as.is=T)
LMISS=read.table("pe9_missing.lmiss", header=T, as.is=T)

png('pe9_missing.png')
oldpar=par(mfrow=c(1,2))
plot( (1:dim(IMISS)[1])/(dim(IMISS)[1]-1), sort(1-IMISS$F_MISS), 
 main="Ordered individual call rate", xlab="Quantile", ylab="Call Rate"); grid()
plot( (1:dim(LMISS)[1])/(dim(LMISS)[1]-1), sort(1-LMISS$F_MISS), 
 main="Ordered SNP coverage", xlab="Quantile", ylab="Coverage" ); grid()
par(oldpar)
dev.off()
