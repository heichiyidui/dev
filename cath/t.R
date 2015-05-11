library(bclust)

cx <- as.matrix(read.table('t2.in'))

mc.cx <- meancss(cx)

optimfunc<-function(theta)
{
-loglikelihood(x.mean=mc.cx$mean,x.css=mc.cx$css,
repno=mc.cx$repno,transformed.par=theta)#compute - log likelihood
}

transpar<-optim(rep(0,6),optimfunc,method="BFGS")$par

bhc=bclust(cx,transformed.par=transpar)

tree_cut=cutree(tree=bhc,h=bhc$cut)

write.table(file='t.out',tree_cut,col.names = FALSE,quote = FALSE)


# ulimit -s unlimited

x   <- as.matrix(read.table('fxf.in'))
clu <- as.matrix(read.table('t.clu'))

kc <- kmeans(x,clu,iter.max = 100)

