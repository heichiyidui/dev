#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

################################################################################
# train network to predict DSSP 3->7 states definitions                        #
################################################################################

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

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
                'on DSSP 3->7 states prediction',\
    epilog=' [-hu number_of_hidden_units -rs random_seed]')

parser.add_argument('-hu', dest='hidden_number',default=20,type=int,\
                    help='number of hidden units, default 20')

parser.add_argument('-rs', dest='random_seed',default=0,type=int,\
                    help='random seed, default 0')

args=parser.parse_args()
random.seed(514+args.random_seed)

################################################################################
# read the data set                                                            #
################################################################################

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

#######################################
# get the DSSP 7 states definitions 

MAP_8_TO_3 = {'X':'X',
              'E':'E','B':'E',
              'H':'H','G':'H',
              'C':'C','S':'C','T':'C','I':'C'}

KEEP_RATE = {'C':0.09439,'D':0.97419,'E':0.27341,'F':0.97422,'G':0.99937,\
             'H':0.13378,'I':1.00000}

MAP_7_TO_I = {'C':0,'D':1,'E':2,'F':3,'G':4,'H':5,'I':6}

states  = {}
str_file=open('../cath/index/cath_s35.dssp') 
for line in str_file:
    dom_id = line.strip()[1:]
    ss_line = str_file.readline().strip()
    if dom_id not in dom_ids: continue
    
    ss_3='' 
    for c in ss_line:
        ss_3 += MAP_8_TO_3[c]

    ss_7=''
    ss_7 += ss_3[0]
    for i in range(len(ss_3)-2):
        ss_char = ss_3[i+1]
        if ss_3[i:i+3] in ('CHH','EHH'):
            ss_char='G'
        if ss_3[i:i+3] in ('HHC','HHE'):
            ss_char='I'
        if ss_3[i:i+3] in ('HEE','CEE'):
            ss_char='D'
        if ss_3[i:i+3] in ('EEH','EEC'):
            ss_char='F'
        ss_7 += ss_char
    ss_7 += ss_3[-1]

    ss=[]
    for c in ss_7:
        if c == 'X': 
            ss.append('X')
            continue
        
        if random.random() > KEEP_RATE[c]:
            ss.append('X')
            continue 
        target = [0.0] *7
        target[MAP_7_TO_I[c]] = 1.0
        ss.append(target)

    states[dom_id]=ss 

str_file.close()

#######################################
# assign training order 
train_res_sum=0
train_index={}
for dom_id in dom_ids:
    t_index=[]
    for i in range(len(states[dom_id])):
        if states[dom_id][i] == 'X': 
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
# training and cross validation                                                #
################################################################################

net=Snn(IN_SIZE,args.hidden_number,7,eta=0.0002)
cv_errors=[999.0]

# training and cross validation iterations
for iter in range(1000):

    # training 
    for dom_id in dom_ids:
        for i in train_index[dom_id]:
            net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
            net.back_propagate(states[dom_id][i])
    
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
                    net.back_propagate(states[dom_id][i])

        # k-test 
        for dom_id in cv_dom_id_lists[k_test]:
            for i in train_index[dom_id]:
                net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
                cv_error += net.ce_error(states[dom_id][i])

    cv_error /= train_res_sum
    cv_errors.append(cv_error)
    net=copy.deepcopy(net_bak)
    print(iter,cv_error)
    if cv_error > cv_errors[-2]:
        break;

net.save_net('y_7_bal_'+str(args.hidden_number)+'_'\
                       +str(args.random_seed)+'.net')
