#!/usr/bin/env python3

#######################################
# read the training and testing set

import random
random.seed()
train_set=[]
test_set=[]

ifile=open('t.in')
for line in ifile:
    cols=line.split()
    if cols[0]=="MCI":
        continue
    
    ab42 = ( float(cols[1]) - 386 ) / 141
    pTau = ( float(cols[3]) - 37.9) / 19.4
    tTau = ( float(cols[5]) - 115 ) / 66.5
    
    input=[ab42, pTau, tTau]
    
    if cols[0]=="AD":
        target=  1.0 
    else:
        target= -1.0 
    
    train_set.append([input,target])

ifile.close()

random.shuffle(train_set)

#######################################
# train the network 

from bpnn import Bpnn
net=Bpnn(3,1)
for i in range(150):
    random.shuffle(train_set)
    train_err=0;
    for j in range(len(train_set)):
        net.propagate     (train_set[j][0])
        net.back_propagate(train_set[j][1])
        train_err += net.se_error(train_set[j][1])
    train_err /= len(train_set)
    print(i,train_err)
    

for pat in train_set:
    print(pat,net.propagate(pat[0]))
