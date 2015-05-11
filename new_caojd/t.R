args = commandArgs(trailingOnly = TRUE)

dom_id=args[1]
dom_len=strtoi(args[2])

c_mat=as.matrix(read.table('t.in',header=FALSE))

png(paste(dom_id,'png',sep='.'),width=dom_len*3+80,height=dom_len*3+130)

image(c_mat, col=(c('#FFFFFF','#000000')), useRaster = TRUE ,axes = FALSE)

axis(1, labels=FALSE, tick=TRUE, at = seq(0, 1,100/dom_len))
axis(2, labels=FALSE, tick=TRUE, at = seq(0, 1,100/dom_len))

axis(1, labels=FALSE, tick=TRUE, at = c(0, 1) )
axis(2, labels=FALSE, tick=TRUE, at = c(0, 1) )
axis(3, labels=FALSE, tick=TRUE, at = c(0, 1) )
axis(4, labels=FALSE, tick=TRUE, at = c(0, 1) )

title(main = dom_id)

dev.off()
