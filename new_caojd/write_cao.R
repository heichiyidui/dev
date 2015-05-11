Pi=scan('cao/Pi')
F=diag(c(t(Pi)))

Q=read.table('cao/Q',header=FALSE)
ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

p_background = array(0,dim=c(400,400))
for (i in 1:400)
    for (j in 1:400)
        p_background[i,j] = Pi[i] * Pi[j]

################################################################################
t=100
p_t = F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)
s_t = log(p_t) - log(p_background)

write.table(s_t,'cao100',col.names=FALSE,row.names=FALSE)

t=150
p_t = F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)
s_t = log(p_t) - log(p_background)

write.table(s_t,'cao150',col.names=FALSE,row.names=FALSE)

t=200 
p_t = F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)
s_t = log(p_t) - log(p_background)

write.table(s_t,'cao200',col.names=FALSE,row.names=FALSE)

t=250 
p_t = F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)
s_t = log(p_t) - log(p_background)

write.table(s_t,'cao250',col.names=FALSE,row.names=FALSE)

################################################################################
# to plot the VTML matrices 

png('cao_250.png',width=1800,height=1720)

par(mar=c(5,4.5,4,7))

image(s_t, axes=FALSE, col=heat.colors(12) )

Labels=c("A","R","N","D","C","Q","E","G","H","I",
         "L","K","M","F","P","S","T","W","Y","V")

axis(side=1, at=(0:19)/20+0.025, labels=Labels)
axis(side=2, at=(0:19)/20+0.025, labels=Labels, las= HORIZONTAL<-1)

library(fields) # for image.plot
image.plot(s_t, col=heat.colors(12), legend.only=T)

dev.off()
