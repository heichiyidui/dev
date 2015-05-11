#! /usr/bin/env python3

# to assign ages back to the t.in table
dfile=open('age.assign')
age={}
for line in dfile:
    age[line.split()[0]]=float(line.split()[1])
dfile.close()

ifile=open('t.in')
for line in ifile:
    id=line.split()[0]
    if id in age.keys():
        print(id,line.split()[1],age[id])
    else:
        print(line[:-1])
ifile.close()
