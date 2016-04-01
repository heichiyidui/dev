#!/usr/bin/Rscript

dom_ids = read.table('t.ls',header=FALSE,as.is=1)

for (dom_id in dom_ids$V1){
    mat_file_name = paste('cont_map/',dom_id,sep='')
    mat = read.table(mat_file_name, header=FALSE)
    mat = as.matrix(mat)

    png_file_name = paste('cont_map_png/',dom_id,'.png',sep='')
    png(png_file_name)

    heatmap(-mat, Rowv = NA, Colv = NA, labRow=NA, labCol=NA,
            margins = c(0.5, 0.5) , scale='none' ,col = topo.colors(5) )

    dev.off()
}
