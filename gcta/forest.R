library(rmeta)
m  <- c(0.115051,0.115309,0.115342,0.115074,0.11783,0.126355,0.116293,0.136268)
se <- c(0.047003,0.030072,0.013719,0.057717,0.015229,0.043272,0.00864,0.015423)
meta.summaries(m, se)

# Summary effect=0.119   95% CI (0.108, 0.131)


m <- c(0.115051,NA,0.115309,NA,0.115342,NA,0.115074,NA,0.11783,NA,0.126355,NA,0.116293,NA,0.136268,NA,NA,0.119,NA,NA)
l <- c(0.022925,NA,0.056368,NA,0.088453,NA,0.001949,NA,0.08798,NA,0.041542,NA,0.099359,NA,0.106039,NA,NA,0.108,NA,NA)
u <- c(0.207177,NA,0.174250,NA,0.142231,NA,0.228199,NA,0.14768,NA,0.211168,NA,0.133227,NA,0.166497,NA,NA,0.131,NA,NA)
tabletext <- cbind(c('NIH',NA,'Ultrech-1',NA,'Ultrech-2',NA,'Ireland',NA,'MGH',NA,'NIH-IT',NA,'SLAGEN',NA,'UK',NA,NA,'Summary \n(Random effects)',NA,NA))
pdf('t1.pdf')
forestplot(tabletext,m,l,u,is.summary=c(rep(FALSE,17),TRUE,FALSE,FALSE),,col=meta.colors(line="darkblue",summary="royalblue"),xticks=c(0,0.1,0.2,0.3),xlab='Heritability, GRM cutoff 0.025')
dev.off()

m <-c(0.1150,NA,0.1150,NA,0.1150,NA,0.1150,NA,0.1187,NA,0.1266,NA,0.1152,NA,0.1415,NA,NA,0.119,NA,NA)
l <-c(0.0292,NA,0.0618,NA,0.0981,NA,0.0021,NA,0.1005,NA,0.0451,NA,0.1014,NA,0.1143,NA,NA,0.110,NA,NA)
u <-c(0.2008,NA,0.1682,NA,0.1319,NA,0.2280,NA,0.1369,NA,0.2081,NA,0.1290,NA,0.1687,NA,NA,0.127,NA,NA)
tabletext <- cbind(c('NIH',NA,'Ultrech-1',NA,'Ultrech-2',NA,'Ireland',NA,'MGH',NA,'NIH-IT',NA,'SLAGEN',NA,'UK',NA,NA,'Summary \n(Random effects)',NA,NA))
pdf('t2.pdf')
forestplot(tabletext,m,l,u,is.summary=c(rep(FALSE,17),TRUE,FALSE,FALSE),,col=meta.colors(line="darkblue",summary="royalblue"),xticks=c(0,0.1,0.2,0.3),xlab='Heritability')
dev.off()
