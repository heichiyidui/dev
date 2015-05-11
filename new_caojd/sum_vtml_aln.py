#!/usr/bin/env python3 
##!/home/klinbrc/bin/python3 

dom_ids=open('index/cath_dom_rep.ls').read().split()

sum_aln_mats=[]
for i in range(100):
    sum_aln_mats.append([0]*400)

ifile=open('t_dis_split.ls')
dis_split_list=[]
for line in ifile:
    dis_split_list.append(float(line.strip()))
ifile.close() # split alignments into 100 different sum matrices 

for dom_id in dom_ids:
    dis_file=open('vtml_dis/'+dom_id)
    aln_file=open('vtml_aln/'+dom_id)
    
    for line in dis_file:
        t_dis = float(line.strip())
        aln_line = aln_file.readline().strip()
        
        if t_dis < 1.5: 
            continue
        
        aln_class=0
        for i in range(len(dis_split_list)):
            if t_dis > dis_split_list[i]:
                aln_class += 1
        
        cols=aln_line.split()
        
        for i in range(400):
            sum_aln_mats[aln_class][i] += int(cols[i])
    
    dis_file.close()
    aln_file.close()

for i in range(100):
    for j in range(400):
        print(sum_aln_mats[i][j],end=' ')
    print()


