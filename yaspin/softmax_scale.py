#!/usr/bin/env python3
import sys
import copy
import pickle
import random
from bpnn import Bpnn
import dssp_to_seven

WING_SIZE = 8*20
WIN_SIZE  = WING_SIZE*2 + 20

################################################################################
# 1. read PSSM and DSSP 

ifile=open('index/test.ls')
test_ids=ifile.read().split()
ifile.close()

pssms=[]
seven_state_defs=[]
targets=[]
residue_lists=[]

for domain_id in test_ids:
    ifile=open('b_in/'+domain_id,'br')
    caln=pickle.load(ifile)
    dssp=pickle.load(ifile)
    acc =pickle.load(ifile)
    pssm=pickle.load(ifile)
    cont_map=pickle.load(ifile)
    ifile.close()
    
    wing_pssm=[-2.0]*WING_SIZE
    pssms.append(wing_pssm + pssm + wing_pssm)
    
    seven_state_def = dssp_to_seven.dssp_to_seven(dssp)
    seven_state_defs.append(seven_state_def)
    
    target=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]*len(dssp)
    residue_list=[]
    for i in range(len(dssp)):
        if seven_state_def[i] == -1:
            continue
        residue_list.append(i)
        target[i*7+seven_state_def[i]]=1.0
        
    targets.append(target)
    residue_lists.append(residue_list)
    
#################################################################################
# 2. test the neural network

domain_list=[]
for i in range(len(test_ids)):
    domain_list.append(i)

net=Bpnn(WIN_SIZE,15,7,eta=1e-6)
net.read_net('t_0.net')

for softmax_scale in range(7,8):    
    
    # to get the average sse error
    avg_ce_err = 0.0;
    sum_res = 0
    for i in domain_list:
        sum_res += len(residue_lists[i])
        for j in residue_lists[i]:
            sum_res+=1
            net.propagate(pssms[i][j*20:j*20 + WIN_SIZE])
            net.softmax(softmax_scale)
            print(*targets[i][j*7:j*7+7],end=' ')
            print(*net.act_o)
            avg_ce_err+=net.ce_error(targets[i][j*7:j*7+7])
    avg_ce_err/=sum_res
    
#    print(softmax_scale,avg_ce_err)


