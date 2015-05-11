
dis = scan('t.dis',what = double(0))
Pi  = scan('Pi',what = double(0))

n_mats = array(scan('t.mat',what = double(0)))

n_mat = array(n_mats, dim = c(400,400,100) )

for (i in 2:100) {
    n_mat[,,i] = n_mat[,,i] + 1
    n_mat[,,i] = n_mat[,,i] / sum(n_mat[,,i])
}

p_mat = array(dim=c(400,400,100))

for (i in 1:100){
    p1=n_mat[,,i]
    p1= diag(1/Pi)  %*% p1
    p_mat[,,i]= p1
}

# normalize the P(t) matrices
for (l in 1:100){
    a= p_mat[,,l]
    b = rowSums(a)
    for (i in 1:400){
        for (j in 1:400){
            p_mat[i,j,l] = p_mat[i,j,l]/b[i]
        }
    }
}

################################################################################
# to get Q
alpha = 0.02 
# R = array(0,dim = c(400,400))
R = exp(-alpha * 0) * p_mat[,,1] * dis[2] * 0.5 

for (i in 2:99){
    t = dis[i]
    t_range = (dis[i+1] - dis[i-1]) * 0.5
    R = R + exp(-alpha * t) * p_mat[,,i] * t_range
}

t = dis[100]
t_range = dis[100] - dis[99]

R = R + exp(-alpha * t) * p_mat[,,100] * t_range 

inv_R = solve(R)

I = diag(400)

Q = alpha * I - inv_R
# sum(Q[] < 0) = 78631

################################################################################
# to normalize Q so that I + Q = PAM1 
# sum(abs(Q-t(Q))) = 3.417407
# sum(abs(Q)) = 11.98895

PiQ= diag(Pi) %*% Q

# sum(abs(PiQ)) = 0.02971705 
# sum(abs(PiQ-t(PiQ))) = 2.523508e-17

PiQ = ( PiQ + t(PiQ) ) * 0.5

PiQ = PiQ - diag(diag(PiQ))

if (min(PiQ) < 0) {
    for (i in 1:400)
        for (j in 1:400)
            if (PiQ[i,j] < 0)
                PiQ[i,j]=0
} # PiQ _ij when i != j should always be non-negative

c1 = 0.01 / sum(PiQ)

PiQ = PiQ * c1

Q2 = diag(1/Pi) %*% PiQ

Q2 = Q2 - diag(diag(Q2)) 

Q2 = Q2 - diag(rowSums(Q2))

write.table(Q2,'Q2',col.names=FALSE,row.names=FALSE)

################################################################################
# to get the new Pi
ev <- eigen(Q2)
Lambda=ev$values
S=ev$vectors
inv_S=solve(S)

p100 = diag(c(t(Pi)))  %*%  (S %*% diag(exp(Lambda *100 )) %*% inv_S)
Pi2=colSums(p100)
write.table(Pi2,'Pi2',col.names=FALSE,row.names=FALSE)

################################################################################
# to plot S matrices
# from http://www.phaget4.org/R/image_matrix.html

# source("http://www.phaget4.org/R/myImagePlot.R") 

# p_background = array(0,dim=c(400,400))
# for (i in 1:400) for (j in 1:400) p_background[i,j] = Pi[i]*Pi[j]
#
# s_100 = log(p100) - log(p_background)
# myImagePlot(s_100)

