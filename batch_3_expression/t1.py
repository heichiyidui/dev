#!/usr/bin/env python3

#######################################
# read the training and testing set

import random
random.seed(514)
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
    
    input=(ab42, pTau, tTau)
    
    if cols[0]=="AD":
        target=(1.0, 0.0)
    else:
        target=(0.0, 1.0)
    
    if random.random()>0.2:
        train_set.append([input,target])
    else:
        test_set.append([input,target])
ifile.close()

random.shuffle(train_set)
random.shuffle(test_set)

#######################################
# train the network 

from snn import Snn
net=Snn(3,2,2,0.01)
for i in range(500):
    ce_error=0
    for j in range(len(train_set)):
        net.propagate(train_set[j][0])
        net.back_propagate(train_set[j][1])
        ce_error += net.ce_error(train_set[j][1])
    print(i,ce_error/len(train_set),end=' ')
    
    ce_error=0
    for j in range(len(test_set)):
        net.propagate(test_set[j][0])
        ce_error += net.ce_error(test_set[j][1])
    print("test:",ce_error/len(test_set))
net.save_net('t.net')