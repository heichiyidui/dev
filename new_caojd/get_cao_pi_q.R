################################################################################
# read matrices

dis = as.matrix(read.table('cao_sum.dis',header=FALSE))

n_mats = array(scan('cao_aln.in',what = double(0)))
n_mat = array(n_mats, dim = c(400,400,100) )

# normalize the P(t) matrices
p_mat = array(dim=c(400,400,100))
for (i in 1:100){
    p1=n_mat[,,i]
    pi1=colSums(p1)
    p1= diag(c(1/pi1))  %*% p1
    p_mat[,,i]= p1
}

################################################################################
# estimate Q

alpha = 0.02 
R = exp(-alpha * dis[1] * 0.25) * diag(400) * dis[1] *0.5

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

Q = alpha * diag(400) - inv_R

################################################################################
# to normalize Q so that I + Q = PAM1 
# sum(abs(Q)) = 9.3
# sum(abs(Q-t(Q))) = 2.4

Pi  = as.matrix(read.table('cao/Pi',header=FALSE))
PiQ= diag(c(Pi)) %*% Q

# sum(abs(PiQ)) = 0.023
# sum(abs(PiQ-t(PiQ))) = 0.00071

PiQ = ( PiQ + t(PiQ) ) * 0.5

PiQ = PiQ - diag(diag(PiQ))

if (min(PiQ) < 0) {
    for (i in 1:400)
        for (j in 1:400)
            if (PiQ[i,j] < 0)
                PiQ[i,j]=0
} # PiQ _ij when i != j should always be non-negative
# 54706 negative numbers in the first run

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
inv_S=solve(ev$vectors)

p10000 = S %*% diag(exp(Lambda *10000 )) %*% inv_S
Pi2 = colSums(p10000)/sum(p10000)

write.table(Pi2,'Pi2',col.names=FALSE,row.names=FALSE)
