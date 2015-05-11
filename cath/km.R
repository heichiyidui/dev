args <- commandArgs(trailingOnly = TRUE)
stem <- args[1]

x <- as.matrix(read.table(stem))
kc=kmeans(x,50,iter.max = 100)
write.table(file=paste(stem,".out",sep=""),kc$centers,
    row.names=FALSE,col.names=FALSE)

# Rscript km.R xaa 
