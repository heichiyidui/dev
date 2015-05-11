#!/usr/bin/env python3

ifile=open('t.in')
seqs=[]
fastaBuffer=''
for line in ifile:
    if line.startswith('>'):
        if fastaBuffer !='':
            seqs.append(fastaBuffer)
        fastaBuffer=line
    else:
        fastaBuffer=fastaBuffer+line

ifile.close()

if fastaBuffer !='':
    seqs.append(fastaBuffer)
    fastaBuffer=''

print(len(seqs))
for seq in seqs:
    print(seq)
    
    