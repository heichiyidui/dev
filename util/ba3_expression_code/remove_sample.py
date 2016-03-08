#!/usr/bin/env python3

to_remove_id_ls =  open('to_remove_sample.ls').read().split()


ifile=open('Genome studio/BATCH3 91014/Groupprobes.txt')
for i in range(8):
    print(ifile.readline(),end='')

title_line_cols=ifile.readline()[:-1].split('\t')
to_remove_ls=[]
new_title_line_cols=[]
for i in range(len(title_line_cols)):
    id_found=False
    for sample_id in to_remove_id_ls:
        if title_line_cols[i].startswith(sample_id):
            id_found=True
            break;
    
    if id_found:
        to_remove_ls.append(i)
    else:
        new_title_line_cols.append(title_line_cols[i])

print('\t'.join(new_title_line_cols))

for line in ifile:
    cols=line[:-1].split('\t')
    new_cols=[]
    for i in range(len(cols)):
        if i in to_remove_ls:
            continue
        new_cols.append(cols[i])
    print('\t'.join(new_cols))
ifile.close()
