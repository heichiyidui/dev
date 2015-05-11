
Pi=scan('vtml/Pi')
F=diag(c(t(Pi)))

Q=read.table('vtml/Q',header=FALSE)
ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

p_background = array(0,dim=c(20,20))
for (i in 1:20)
    for (j in 1:20)
        p_background[i,j] = Pi[i] * Pi[j]

################################################################################
t=50 
p_t = F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)
s_t = log(p_t) - log(p_background)

write.table(s_t,'vtml50',col.names=FALSE,row.names=FALSE)

t=75 
p_t = F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)
s_t = log(p_t) - log(p_background)

write.table(s_t,'vtml75',col.names=FALSE,row.names=FALSE)

t=100 
p_t = F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)
s_t = log(p_t) - log(p_background)

write.table(s_t,'vtml100',col.names=FALSE,row.names=FALSE)

t=150 
p_t = F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)
s_t = log(p_t) - log(p_background)

write.table(s_t,'vtml150',col.names=FALSE,row.names=FALSE)

################################################################################
# to plot the VTML matrices 

png('vtml_150.png',width=800,height=720)

par(mar=c(5,4.5,4,7))

image(s_t, axes=FALSE, col=heat.colors(12) )

Labels=c("A","R","N","D","C","Q","E","G","H","I",
         "L","K","M","F","P","S","T","W","Y","V")

axis(side=1, at=(0:19)/19, labels=Labels)
axis(side=2, at=(0:19)/19, labels=Labels, las= HORIZONTAL<-1)

library(fields) # for image.plot
image.plot(s_t, col=heat.colors(12), legend.only=T)

dev.off()
