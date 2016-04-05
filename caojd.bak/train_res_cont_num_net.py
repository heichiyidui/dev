#!/usr/bin/env python3
##!/home/klinbrc/bin/python3
import sys 
while '/share/python/lib64/python' in sys.path :
    sys.path.remove('/share/python/lib64/python')
from numpy import *

################################################################################
# train network to predict 3 states of the global contact number distributions #
################################################################################

CONT_TARGETS=array([
[ 8.476546e-01, 1.502704e-01, 2.075020e-03 ], 
[ 3.718094e-01, 6.050626e-01, 2.312803e-02 ], 
[ 5.708059e-02, 8.526955e-01, 9.022389e-02 ], 
[ 5.608695e-03, 7.691180e-01, 2.252733e-01 ], 
[ 4.385163e-04, 5.520042e-01, 4.475572e-01 ], 
[ 2.667320e-05, 3.082175e-01, 6.917558e-01 ], 
[ 1.307043e-06, 1.386429e-01, 8.613557e-01 ], 
[ 5.643476e-08, 5.495155e-02, 9.450484e-01 ], 
[ 2.301705e-09, 2.057351e-02, 9.794265e-01 ], 
[ 9.178670e-11, 7.531200e-03, 9.924688e-01 ], 
[ 3.629603e-12, 2.733815e-03, 9.972662e-01 ], 
[ 1.430880e-13, 9.893246e-04, 9.990107e-01 ], 
[ 5.634598e-15, 3.576218e-04, 9.996424e-01 ], 
[ 2.217927e-16, 1.292212e-04, 9.998708e-01 ], 
[ 8.729077e-18, 4.668533e-05, 9.999533e-01 ] ])

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

import copy
import random
from snn import Snn

################################################################################
# parse the arguments                                                          #
################################################################################

import argparse

parser = argparse.ArgumentParser(\
    description='6-fold cross-validation training of neural network '+\
                'on 3-states global contact number distribution',\
    epilog=' [-hu number_of_hidden_units -rs random_seed ]')

parser.add_argument('-hu', dest='hidden_number',  default=60,type=int,\
                    help='number of hidden units, default 60')

parser.add_argument('-rs', dest='random_seed',    default=0, type=int,\
                    help='random seed, default 0')

args=parser.parse_args()
random.seed(514+args.random_seed)


################################################################################
# read the data set                                                            #
################################################################################

dom_ids = open('index/s.ls').read().split()
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
# read domain sequences

dom_seqs={}
ifile=open('../cath/index/cath_s35.cseq')
for id_line in ifile:
    dom_id = id_line[1:].strip()
    seq    = ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    dom_seqs[dom_id]=seq 
ifile.close()

dom_num_2_seq_num={}
for dom_id in dom_ids:
    cseq=dom_seqs[dom_id]

    num_2_seq_num={}
    dom_index=0
    for i in range(len(cseq)):
        if cseq[i].isupper():
            num_2_seq_num[dom_index]=i
            dom_index+=1

    dom_num_2_seq_num[dom_id]=num_2_seq_num

#######################################
# read contact maps 

dom_targets={}
ifile=open('index/cath_s35.condef')
for line in ifile:
    c_line=ifile.readline()
    dom_id = line.split()[0][1:]
    if dom_id not in dom_ids: continue 
    
    dom_len = int (line.split()[1])
    cont_num = int (line.split()[2])
    cols=c_line.split()
    
    cont_map=zeros((dom_len,dom_len))
    
    for i in range(cont_num):
        c_i=int(cols[i*2])
        c_j=int(cols[i*2+1])
        cont_map[c_i,c_j]=1
        cont_map[c_j,c_i]=1
    
    targets=['NA'] * len(dom_seqs[dom_id])
    for i in range(dom_len):
        seq_index=dom_num_2_seq_num[dom_id][i]
        targets[seq_index]= CONT_TARGETS[int(sum(cont_map[i]))]
    
    dom_targets[dom_id]=targets 
    
ifile.close()
    
#######################################
# training index     
train_index={}
train_res_sum=0
for dom_id in dom_ids:
    t_index=[]
    for i in range(len(dom_targets[dom_id])):
        if dom_targets[dom_id][i]=='NA':
            continue 
        t_index.append(i)
        
    random.shuffle(t_index)
    train_index[dom_id]=t_index 
    train_res_sum += len(t_index)

#######################################
# 6 fold cross validation, domain index lists 
cv_dom_id_lists=[[],[],[], [],[],[]]
for dom_id in dom_ids:
    cv_dom_id_lists[random.randint(0,5)].append(dom_id)


################################################################################
# training and cross validation                                                #
################################################################################
net=Snn(IN_SIZE,args.hidden_number,3,eta=0.0002)
cv_errors=[999.0]

# training and cross validation iterations
for iter in range(200):

    # training 
    for dom_id in dom_ids:
        for i in train_index[dom_id]:
            net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
            net.back_propagate(dom_targets[dom_id][i])
    
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
                    net.back_propagate(dom_targets[dom_id][i])

        # k-test 
        for dom_id in cv_dom_id_lists[k_test]:
            for i in train_index[dom_id]:
                net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
                cv_error += net.ce_error(dom_targets[dom_id][i])

    cv_error /= train_res_sum
    cv_errors.append(cv_error)
    net=copy.deepcopy(net_bak)
    print(iter,cv_error)
    if cv_error > cv_errors[-2]:
        break;

net.save_net('c_3_'+str(net.nh-1)+'_'+'{:02d}'.format(args.random_seed)\
             +'.net')
################################################################################
# the end                                                                      #
################################################################################

