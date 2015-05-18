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
        target=[1.0 ]
    else:
        target=[-1.0 ]
    
    if random.random()>0.2:
        train_set.append([input,target])
    else:
        test_set.append([input,target])
ifile.close()

random.shuffle(train_set)
random.shuffle(test_set)

#######################################
# train the network 

from bpnn import NN
net=NN(3,1,1)
for i in range(3):
    random.shuffle(train_set)
    net.train(train_set,iterations=1,N=0.001)
    train_err=0;
    for pat in train_set:
        prediction = net.update(pat[0])
        train_err+= 0.5*(prediction[0]-pat[1][0])**2
    train_err /= len(train_set)
    print("train_err:",train_err)

test_err=0
for pat in test_set:
    prediction = net.update(pat[0])
    test_err += 0.5 * (prediction[0] - pat[1][0]) ** 2
test_err /=len(test_set)
print("test_err:",test_err)
net.test(test_set)
    
#     ce_error=0
#     for j in range(len(test_set)):
#         net.propagate(test_set[j][0])
#         ce_error += net.ce_error(test_set[j][1])
#     print("test:",ce_error/len(test_set))
# net.save_net('t.net')