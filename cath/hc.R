args <- commandArgs(trailingOnly = TRUE)
stem <- args[1]

cx <- as.matrix(read.table(stem))
d <- dist(cx, method = "manhattan") # distance matrix
fit <- hclust(d, method="ward")
groups <- cutree(fit, k= 80)
write.table(file=paste(stem,".out",sep=""),groups,
    row.names=FALSE,col.names=FALSE)

# Rscript hc.R xaa 
