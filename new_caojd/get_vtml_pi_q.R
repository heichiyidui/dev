################################################################################
# to read the matrices

dis = as.matrix(read.table('vtml_sum.dis',header=FALSE))

n_mats = array(scan('vtml_aln.in',what = double(0)))
n_mat = array(n_mats, dim = c(20,20,100) )

# Pi=array(0,dim=c(20))
# for (i in 1:100){
#     Pi=Pi+colSums(n_mat[,,i])
# }
# Pi = Pi / sum(Pi)

# normalize the P(t) matrices
p_mat = array(dim=c(20,20,100))
for (i in 1:100){
    p1=n_mat[,,i]
    pi1=colSums(p1)
    p1= diag(c(1/pi1))  %*% p1
    p_mat[,,i]= p1
}

################################################################################
# to get Q
alpha = 0.020 
R = exp(-alpha * dis[1] * 0.25) * diag(20) * dis[1] *0.5

R = R + exp(-alpha * dis[1]) * p_mat[,,1] * dis[1] 

for (i in 2:99){
    t = dis[i]
    t_range = (dis[i+1] - dis[i-1]) * 0.5
    R = R + exp(-alpha * t) * p_mat[,,i] * t_range
}

t = dis[100]
t_range = dis[100] - dis[99]

R = R + exp(-alpha * t) * p_mat[,,100] * t_range 

inv_R = solve(R)

Q = alpha * diag(20) - inv_R

################################################################################
# to normalize Q so that I + Q = PAM1 
# sum(abs(Q-t(Q))) = 0.095
# sum(abs(Q)) = 0.43

Pi  = as.matrix(read.table('vtml/Pi',header=FALSE))

PiQ= diag(c(Pi)) %*% Q

# sum(abs(PiQ))= 0.021
# sum(abs(PiQ-t(PiQ))) = 0.00097

PiQ = ( PiQ + t(PiQ) ) * 0.5

PiQ = PiQ - diag(diag(PiQ))

if (min(PiQ) < 0) {
    for (i in 1:20)
        for (j in 1:20)
            if (PiQ[i,j] < 0)
                PiQ[i,j]=0
} # PiQ _ij when i != j should always be non-negative

PiQ = PiQ * (0.01 / sum(PiQ))

Q2 = diag(1/c(Pi)) %*% PiQ

Q2 = Q2 - diag(rowSums(Q2))

write.table(Q2,'Q2',col.names=FALSE,row.names=FALSE)

################################################################################
# to get the new Pi
Q=Q2
ev <- eigen(Q)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

p10000 = S %*% diag(exp(Lambda *10000 )) %*% inv_S
Pi2 = colSums(p10000)/sum(p10000)

write.table(Pi2,'Pi2',col.names=FALSE,row.names=FALSE)
