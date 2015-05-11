library(pracma,lib='~/R')
#library(pracma)

Pi=as.matrix(read.table('cao/Pi',header=FALSE))
F=diag(c(Pi))
Q=as.matrix(read.table('cao/Q',header=FALSE))

ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

args = commandArgs(trailingOnly = TRUE)
domain_ids=scan(args[1],'')
# domain_ids=scan('t.ls','')

for (l in 1:length(domain_ids) ){
    domain_id = domain_ids[l]
    
    n_mats = scan(paste('cao_aln/',domain_id,sep='/'),what = double(0))
    
    num_aln = length(n_mats) / 160000
    list_dis=c()
    
    for (i in 1:num_aln){
        j=i*160000
        n_mat = array(n_mats[(i*160000-159999):(i*160000)],dim=c(400,400))
        
        loglikeN <- function (t) {
             -sum(n_mat * log(F %*% (S %*% diag(exp(Lambda *t )) %*% inv_S)) )
        }
        
        min_dis=fminbnd(loglikeN,1,1000)
        list_dis=c(list_dis,min_dis$xmin)
        
    }
    
    write.table(list_dis,file = paste('cao_dis',domain_id,sep='/'),
                row.names = FALSE, col.names = FALSE)
#    break;
}
