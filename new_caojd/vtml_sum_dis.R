#library(pracma,lib='~/R')
library(pracma)

Pi=as.matrix(read.table('vtml/Pi',header=FALSE))
F=diag(c(Pi))

Q=as.matrix(read.table('vtml/Q',header=FALSE))
ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

n_mats = as.matrix(read.table('vtml_aln.in',header=FALSE))
list_dis=c()

for (i in 1:100){
    n_mat = array(n_mats[i,], dim=c(20,20))
    
    # using local Pi diagnal matrix makes almost no difference
    # Pi2=as.matrix(colSums(n_mat)/sum(n_mat))
    # F=diag(c(Pi2)) 
    
    loglikeN <- function (t) {
         -sum(n_mat * log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)) )
    }
    
    min_dis=fminbnd(loglikeN,1,1000)
    list_dis=c(list_dis,min_dis$xmin)
}
write.table(list_dis,file='vtml_sum.dis',row.names = FALSE, col.names = FALSE)

