#library(pracma,lib='~/R')
library(pracma)

Pi=as.matrix(read.table('vtml/Pi',header=FALSE))
F=diag(c(Pi))
Q=as.matrix(read.table('vtml/Q',header=FALSE))

ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

args = commandArgs(trailingOnly = TRUE)
n_mat= array(scan(args[1],what = double(0)), dim=c(20,20))

loglikeN <- function (t) {
    -sum(n_mat * log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)) )
}
        
min_dis=fminbnd(loglikeN,1,1000)
print(min_dis$xmin)