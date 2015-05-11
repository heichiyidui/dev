#! /usr/bin/env python3

ifile=open('t.ls')
probRs=[]
for line in ifile:
    probRs.append(line.strip())
ifile.close()

ifile=open('t.in')
for line in ifile:
    rsid=line.split()[0]
    
    if rsid in probRs:
        v1=1/float(line.split()[1])
        v2=1/float(line.split()[3])
        v3=1/float(line.split()[2])
        print(rsid,v1,v2,v3)
    else:
        print(line[:-1])
ifile.close()
