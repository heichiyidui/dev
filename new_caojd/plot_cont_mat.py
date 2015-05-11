#!/usr/bin/env python3 
import subprocess

dom_ids=open('index/cath_dom_rep.ls').read().split()

dom_lens={}
ifile=open('index/cath_s35.seg')
for line in ifile:
    cols=line.split()
    if cols[0] not in dom_ids:
        continue
    dom_lens[cols[0]]=int(cols[1])
ifile.close()

ifile=open('index/cath_s35.condef')
for line in ifile:
    dom_id = line[1:-1]
    c_line=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    
    dom_len=dom_lens[dom_id]
    cont_mat=[]
    for i in range(dom_len):
        cont_mat.append([0] * dom_len )
    
    cols=c_line.split()
    for col in cols:
        
        (str_i,str_j)=col.split('-')
        i=int(str_i)
        j=int(str_j)
        cont_mat[i][j]=1
        cont_mat[j][i]=1
    
    ofile=open('t.in','w')
    for i in range(dom_len):
        for j in range(dom_len):
            print(cont_mat[i][j],end=' ',file=ofile)
        print(file=ofile)
    ofile.close()
    
    subprocess.call('Rscript t.R '+dom_id+' '+str(dom_len),shell=True)
    
    #break;
    
ifile.close()

# t.R
#--------------------------------------
# args = commandArgs(trailingOnly = TRUE)
# dom_id=args[1]
# dom_len=strtoi(args[2])
#
# c_mat=as.matrix(read.table('t.in',header=FALSE))
# png(paste(dom_id,'png',sep='.'),width=dom_len*3+80,height=dom_len*3+130)
#
# image(c_mat, col=(c('#FFFFFF','#000000')), useRaster = TRUE ,axes = FALSE)
# axis(1, labels=FALSE, tick=TRUE, at = seq(0, 1,100/dom_len))
# axis(2, labels=FALSE, tick=TRUE, at = seq(0, 1,100/dom_len))
#
# axis(1, labels=FALSE, tick=TRUE, at = c(0, 1) )
# axis(2, labels=FALSE, tick=TRUE, at = c(0, 1) )
# axis(3, labels=FALSE, tick=TRUE, at = c(0, 1) )
# axis(4, labels=FALSE, tick=TRUE, at = c(0, 1) )
#
# title(main = dom_id)
# dev.off()
#
#--------------------------------------
# end of t.R 
