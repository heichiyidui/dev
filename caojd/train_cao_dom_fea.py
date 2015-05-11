#!/usr/bin/env python3
################################################################################
# train networks to predict 2-states CAO contact definition                    #
################################################################################

import sys
import copy
import random
from snn import Snn

################################################################################
# parse the arguments                                                          #
################################################################################

import argparse

parser = argparse.ArgumentParser(\
    description='6-fold cross-validation training of neural network '+\
                'on CAO residue contact prediction',\
    epilog=' [-hu number_of_hidden_units -rs random_seed]')

parser.add_argument('-hu', dest='hidden_number',  default=6,type=int,\
                    help='number of hidden units, default 6')

parser.add_argument('-rs', dest='random_seed',    default=0, type=int,\
                    help='random seed, default 0')

args=parser.parse_args()
random.seed(514+args.random_seed)

################################################################################
# read inputs and targets                                                      #
################################################################################

inputs=[]
targets=[]

ifile=open('train.in')
for line in ifile:
    cols=line.split()
    inputs.append([float(cols[0]),float(cols[1]),float(cols[2]),float(cols[3])])
    if cols[4]=='1':
        targets.append([1.0,0.0])
    else:
        targets.append([0.0,1.0])
    
ifile.close()

# 6 fold cross validation, domain index lists 
train_size=len(inputs)
cv_lists=[[],[],[], [],[],[]]
for i in range(train_size):
    cv_lists[random.randint(0,5)].append(i)

################################################################################
# training and cross validation                                                #
################################################################################

net=Snn(4,args.hidden_number,2,eta=0.0002)

cv_errors=[999.0]
for iter in range(1000):
    # training 
    for i in range(train_size):
        net.propagate(inputs[i])
        net.back_propagate(targets[i])
    
    if iter < 10: # first 10 iterations are just warming up 
        print(iter)
        continue 
    
    if iter %3 > 0 : # cross validation every 3 iterations of training
        print(iter)
        continue 
    
    # cross validation 
    net_bak=copy.deepcopy(net)
    cv_error=0.0
    for k_test in range(6):
        net=copy.deepcopy(net_bak)
        # k-train 
        for k_train in range(6):
            if k_train == k_test: continue

            for i in cv_lists[k_train]:
                net.propagate(inputs[i])
                net.back_propagate(targets[i])

        # k-test 
        for i in cv_lists[k_test]:
            net.propagate(inputs[i])
            cv_error += net.ce_error(targets[i])

    cv_error /= train_size
    cv_errors.append(cv_error)
    net=copy.deepcopy(net_bak)
    print(iter,cv_error)
    if cv_error > cv_errors[-2]:
        break;

net.save_net('c_2_'+str(net.nh-1)+'_'+str(args.random_seed)+'.net')

################################################################################
# end                                                                          #
################################################################################
