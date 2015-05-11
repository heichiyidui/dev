#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

################################################################################
# train network to predict 10 states definitions                               #
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
                'on DSSP 10 states prediction',\
    epilog=' [-hu number_of_hidden_units -rs random_seed -net in_net_filename]')

parser.add_argument('-hu', dest='hidden_number',  default=20,type=int,\
                    help='number of hidden units, default 20')

parser.add_argument('-rs', dest='random_seed',    default=0, type=int,\
                    help='random seed, default 0')

args=parser.parse_args()
random.seed(514+args.random_seed)

################################################################################
# read the data set                                                            #
################################################################################

dom_ids = open('index/train.ls').read().split()
while '' in dom_ids: dom_ids.remove('')
random.shuffle(dom_ids)

#######################################
# read pssm 
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
# read residue exposure
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
        else :
            expo.append(float(exp_str))
    
    expos[dom_id]=expo 

exp_file.close()

#######################################
# read secondary structure 
ss_states  = {}

MAP_8_TO_4={'X':'X',
            'E':'E','B':'E',
            'H':'H',
            'G':'G',
            'C':'C','S':'C','T':'C','I':'C'}

str_file=open('../cath/index/cath_s35.dssp') 
for line in str_file:
    dom_id = line.strip()[1:]
    ss_line = str_file.readline().strip()
    if dom_id not in dom_ids: continue
    
    ss=[]
    for c in ss_line:
        ss.append(MAP_8_TO_4[c])
    
    ss_5=[]
    ss_5.append(ss[0])
    for i in range(0,len(ss)-2):
        if ss[i]+ss[i+1]+ss[i+2] == 'CHH': ss_5.append('Hb'); continue
        if ss[i]+ss[i+1]+ss[i+2] == 'EHH': ss_5.append('Hb'); continue

        ss_5.append(ss[i+1])

    ss_5.append(ss[-1])
    ss_states[dom_id]=ss_5
    
str_file.close()

#######################################
# get traning targets
targets={}

MAP_5_TO_i={'C':0,'E':1,'Hb':2,'H':3,'G':4}

for dom_id in dom_ids:
    states=[]
    for i in range(len(ss_states[dom_id])):
        ss  = ss_states[dom_id][i]
        expo= expos[dom_id][i]
        if ss == 'X' or expo == 'NA':
            states.append('X')
            continue
        state=MAP_5_TO_i[ss]
        if expo < 0.5:
            state +=5
        states.append(state)
    
    targets[dom_id]=[]
    for state in states:
        if state=='X':
            targets[dom_id].append('X')
            continue
        target = [0.0] * 10
        target[state]=1.0

        targets[dom_id].append(target)


train_res_sum=0
train_index={}
for dom_id in dom_ids:
    t_index=[]
    for i in range(len(targets[dom_id])):
        if targets[dom_id][i] == 'X': 
            continue
        t_index.append(i)
    
    random.shuffle(t_index)
    train_index[dom_id]=t_index
    train_res_sum += len(t_index)

# 6 fold cross validation, domain index lists 
cv_dom_id_lists=[[],[],[], [],[],[]]
for dom_id in dom_ids:
    cv_dom_id_lists[random.randint(0,5)].append(dom_id)

# ################################################################################
# # training and cross validation                                                #
# ################################################################################

net=Snn(IN_SIZE,args.hidden_number,10,eta=0.0002)

cv_errors=[999.0]

# training and cross validation iterations
for iter in range(1000):

    # training 
    for dom_id in dom_ids:
        for i in train_index[dom_id]:
            net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
            net.back_propagate(targets[dom_id][i])
    
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
                    net.back_propagate(targets[dom_id][i])

        # k-test 
        for dom_id in cv_dom_id_lists[k_test]:
            for i in train_index[dom_id]:
                net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
                cv_error += net.ce_error(targets[dom_id][i])

    cv_error /= train_res_sum
    cv_errors.append(cv_error)
    net=copy.deepcopy(net_bak)
    print(iter,cv_error)
    if cv_error > cv_errors[-2]:
        break;

net.save_net('y_10_'+str(net.nh-1)+'_'+'{:02d}'.format(args.random_seed)+'.net')
