library('bclust',lib.loc='~/R')

cx <- as.matrix(read.table('txt.cents'))

mc.cx <-meancss(cx)

optimfunc<-function(theta)
{
-loglikelihood(x.mean=mc.cx$mean,x.css=mc.cx$css,
repno=mc.cx$repno,transformed.par=theta)#compute - log likelihood
}

transpar<-optim(rep(0,6),optimfunc,method="BFGS")$par

bhc=bclust(cx,transformed.par=transpar)

tree_cut=cutree(tree=bhc,h=bhc$cut)

write.table(file='txt.cents',tree_cut,col.names = FALSE,quote = FALSE)


# ulimit -s unlimited
# /share/apps/R_3.0.2/bin/R
# plot(bhc)

# clues

library(clues)

cx <- as.matrix(read.table('sxs.cents'))
c=clues(cx,strengthMethod='CH')
write.table(file='sxs_clues.clu',c$mem,row.names=FALSE,col.names=FALSE)

cx <- as.matrix(read.table('nxn.cents'))
d=clues(cx,strengthMethod='CH')
write.table(file='nxn_clues.clu',d$mem,row.names=FALSE,col.names=FALSE)
