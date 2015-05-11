#! /usr/bin/env python3

ifile=open('pe3_part.fam')
fam={} # family ids, without the leading 'F', as numbers
for line in ifile:
    fam[line.split()[1]]=int(line.split()[0][1:])
ifile.close()

for i in range(10):
    ifile=open('pe3_part.fam')
    for line in ifile:
        col=line.split()
        id=col[1]
        father=col[2]; mother=col[3]

        if father != '0':
            if fam[father] < fam[id]:
                fam[id]=fam[father]
            if fam[id] < fam[father]:
                fam[father]=fam[id]
        if mother !='0':
            if fam[mother] < fam[id]:
                fam[id]=fam[mother]
            if fam[id] < fam[mother]:
                fam[mother]=fam[id]
    ifile.close()

ifile=open('pe3_part.fam')
for line in ifile:
    col=line.split()
    id=col[1];
    print('F'+str(fam[id]),id,col[2],col[3],col[4],col[5])
ifile.close()
