
dis = scan('t.dis',what = double(0))
Pi  = scan('Pi',what = double(0))
n_mats = array(scan('t.mat',what = double(0)))

n_mat = array(n_mats, dim = c(20,20,100) )

p_mat = array(dim=c(20,20,100))

for (i in 1:100){
    p1=n_mat[,,i]
    p1= diag(1/Pi)  %*% p1
    p_mat[,,i]= p1
}

# normalize the P(t) matrices
for (l in 1:100){
    a= p_mat[,,l]
    b = rowSums(a)
    for (i in 1:20){
        for (j in 1:20){
            p_mat[i,j,l] = p_mat[i,j,l]/b[i]
        }
    }
}

################################################################################
# to get Q
alpha = 0.02 
# R = array(0,dim = c(20,20))
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

I = diag(20)

Q = alpha * I - inv_R

################################################################################
# to normalize Q so that I + Q = PAM1 
# sum(abs(Q-t(Q))) = 0.08386216
# sum(abs(Q)) = 0.379

PiQ= diag(Pi) %*% Q

# sum(abs(PiQ))= 0.01940659
# sum(abs(PiQ-t(PiQ))) = 0.0003342964

PiQ = ( PiQ + t(PiQ) ) * 0.5

PiQ = PiQ - diag(diag(PiQ))

if (min(PiQ) < 0) {
    for (i in 1:20)
        for (j in 1:20)
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
# ----- Define a function for plotting a matrix ----- #
# from http://www.phaget4.org/R/image_matrix.html

# source("http://www.phaget4.org/R/myImagePlot.R") 

#myImagePlot(Q,
# xLabels=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'),
# yLabels=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
#)

