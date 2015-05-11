
Pi=scan('Pi')
F=diag(c(t(Pi)))

Q=read.table('Q',header=FALSE)
ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

p_background = array(0,dim=c(20,20))
for (i in 1:20)
    for (j in 1:20)
        p_background[i,j] = Pi[i] * Pi[j]

p_background = log(p_background)

t=50 
p_t = log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S))
s_t = p_t - p_background 

write.table(s_t,'vtml50',col.names=FALSE,row.names=FALSE)

t=75 
p_t = log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S))
s_t = p_t - p_background 

write.table(s_t,'vtml75',col.names=FALSE,row.names=FALSE)

t=100 
p_t = log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S))
s_t = p_t - p_background 

write.table(s_t,'vtml100',col.names=FALSE,row.names=FALSE)

