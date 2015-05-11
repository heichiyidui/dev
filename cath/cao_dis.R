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

args = commandArgs()
domain_ids=scan(args[3],'')
# domain_ids=scan('t.ls','')

for (i in 1:length(domain_ids) ){
    domain_id = domain_ids[i]
    n_mats = scan(paste('n_mat',domain_id,sep='/'),what = double(0))
    
    num_aln = length(n_mats) / 400 
    list_dis=c()
    for (i in 1:num_aln){
        j=i*400
        n_mat = array(n_mats[(i*400-399):(i*400)],dim=c(20,20))
        
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

