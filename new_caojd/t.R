library(pracma)

Pi=as.matrix(read.table('vtml/Pi',header=FALSE))
F=diag(c(Pi))
Q=as.matrix(read.table('vtml/Q',header=FALSE))

ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

aln_mat = as.matrix(read.table('t.in'))

# t=60

loglikeN <- function (t) {
    -sum(n_mat * log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)) )
}

min_dis=fminbnd(loglikeN,1,1000)
# min_dis$xmin = 91.97997
