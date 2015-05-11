library(pracma)
# on the cluster 
# library(pracma,lib='~/R')

Pi=read.table('Pi',header=FALSE)
F=diag(c(t(Pi)))

Q=read.table('Q',header=FALSE)
ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

n_mats = scan('t.mat',what = double(0))
num_aln = 100 
list_dis=c(0)
for (i in 2:num_aln){
    j=i*160000
    n_mat = array(n_mats[(j-159999):j],dim=c(400,400))
    n_mat = n_mat + 1
    
    loglikeN <- function (t) {
         -sum(n_mat * log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)) )
    }
    
    min_dis=fminbnd(loglikeN,1,1000,maxiter=50)
    list_dis=c(list_dis,min_dis$xmin)
}
write.table(list_dis,file='t.dis',row.names = FALSE, col.names = FALSE)

