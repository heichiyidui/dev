#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

################################################################################
# test the snn neural network on residue exposure prediction                   #
# 1. test the continue training targets vs discrete training targets.          #
# 2. test the different numbers of hidden units.                               #
################################################################################

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

import sys
import copy
import random
from snn import Snn

################################################################################
# parse the arguments

import argparse

parser = argparse.ArgumentParser(\
    description='6-fold cross-validation training of neural network '+\
                'on residue exposure prediction',\
    epilog=' [-hu number of hidden units] [-bt] [-rs random seed]')

parser.add_argument('-hu', dest='hidden_number',default=10,type=int,\
                    help='number of hidden units, default 10')
parser.add_argument('-rs', dest='random_seed',default=0,type=int,\
                    help='random seed, default 0')

parser.add_argument('-bt', dest='is_binary_target',action='store_true',\
                    help='training targets is binary, default False')
parser.set_defaults(is_binary_target=False)

args=parser.parse_args()

random.seed(514+args.random_seed)

################################################################################
# read the training set 
dom_ids = open('index/train.ls').read().split()
while '' in dom_ids: dom_ids.remove('')

random.shuffle(dom_ids)

pssms={}
pssm_file=open('../cath/index/cath_s35.pssm')
for line in pssm_file:
    dom_id = line.strip()[1:]
    pssm_line = pssm_file.readline().strip()
    
    if dom_id not in dom_ids: continue
    
    pssm  = [-2.0] * 20 * WING_SIZE
    
    for residue_pssm in pssm_line.split():
        for score in residue_pssm.split('|'):
            pssm.append(float(score)) 
    
    pssm += [-2.0] * 20 * WING_SIZE
    pssms[dom_id]=pssm 

pssm_file.close()

expos={}
exp_file = open('index/cath_s35.exp')
for line in exp_file:
    dom_id = line.strip()[1:]
    exp_line = exp_file.readline().strip()
    
    if dom_id not in dom_ids: continue 
    
    expo=[]
    for exp_str in exp_line.split():
        if exp_str == 'NA' : 
            expo.append(exp_str)
        else: 
            if args.is_binary_target:
                expo.append([int(float(exp_str)>0.5),int(float(exp_str)<=0.5)])
            
            else :
                expo.append([float(exp_str),1-float(exp_str)])
    
    
    expos[dom_id]=expo 

exp_file.close()

train_res_sum=0
train_index={}
for dom_id in dom_ids:
    t_index=[]
    for i in range(len(expos[dom_id])):
        if expos[dom_id][i] == 'NA': 
            continue
        t_index.append(i)
    
    random.shuffle(t_index)
    train_index[dom_id]=t_index
    train_res_sum += len(t_index)


# 6 fold cross validation, domain index lists 
cv_dom_id_lists=[[],[],[], [],[],[]]
for dom_id in dom_ids:
    cv_dom_id_lists[random.randint(0,5)].append(dom_id)

################################################################################
# training and cross validation

net=Snn(IN_SIZE,args.hidden_number,2,eta=0.0002)

cv_errors=[999.0]
# training and cross validation iterations
for iter in range(1000):

    # training 
    for dom_id in dom_ids:
        for i in train_index[dom_id]:
            net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
            net.back_propagate(expos[dom_id][i])
    
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

            for dom_id in cv_dom_id_lists[k_train]:
                for i in train_index[dom_id]:
                    net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
                    net.back_propagate(expos[dom_id][i])

        # k-test 
        for dom_id in cv_dom_id_lists[k_test]:
            for i in train_index[dom_id]:
                net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
                cv_error += net.ce_error(expos[dom_id][i])

    cv_error /= train_res_sum
    cv_errors.append(cv_error)
    net=copy.deepcopy(net_bak)
    print(iter,cv_error)
    if cv_error > cv_errors[-2]:
        break;

net.save_net('expos_'+str(args.hidden_number)+'_'+str(args.random_seed)+'_'\
    +str(args.is_binary_target)+'.net')

