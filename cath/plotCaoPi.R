pi_cao = array(scan('Pi',what = double(0)),dim=c(20,20))
pi_cao = t(pi_cao)
pi_cao = log(pi_cao)

Pi=scan('vtml/Pi')
p_background = array(0,dim=c(20,20))
for (i in 1:20)
    for (j in 1:20)
        p_background[i,j] = Pi[i] * Pi[j]
p_background = log(p_background)

pi_cao = pi_cao - p_background

source("http://www.phaget4.org/R/myImagePlot.R") 
myImagePlot(pi_cao,
 xLabels=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'),
 yLabels=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
)

png('t.png',width=1000,height=800)
myImagePlot(pi_cao,title='ln( Pi_CAO / Pi_i x Pi_j )',
 xLabels=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'),
 yLabels=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
)
dev.off()

