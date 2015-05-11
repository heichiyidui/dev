library(bclust)

# to remove the limit of stack size
# ulimit -s unlimited
# head -n 25000 fxf.in > t.in 

cx <- as.matrix(read.table('t.in'))
mc.cx <-meancss(cx)

optimfunc<-function(theta){
-loglikelihood(x.mean=mc.cx$mean,x.css=mc.cx$css,
repno=mc.cx$repno,transformed.par=theta)#compute - log likelihood
}

transpar<-optim(rep(0,6),optimfunc,method="BFGS")$par
bhc=bclust(cx,transformed.par=transpar)
tree_cut=cutree(tree=bhc,h=bhc$cut)
write.table(file='t.out',tree_cut,col.names = FALSE,quote = FALSE)

