#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

################################################################################
# test networks on predicting DSSP 3->7 states definitions                     #
################################################################################

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

import sys
import copy
import random
from snn import Snn

net_file_name=sys.argv[1]

################################################################################
# read the data set                                                            #
################################################################################

dom_ids = open('index/test.ls').read().split()
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

import dssp_to_seven
states  = {}
states_int={}
str_file=open('../cath/index/cath_s35.dssp') 
for line in str_file:
    dom_id = line.strip()[1:]
    ss_line = str_file.readline().strip()
    if dom_id not in dom_ids: continue
    
    ss_states = dssp_to_seven.dssp_to_seven(ss_line)
    states_int[dom_id]=ss_states
    ss=[]
    for state in ss_states:
        if state == -1: 
            ss.append('X')
            continue

        target = [0.0] *7
        target[state]=1.0
        ss.append(target)

    states[dom_id]=ss 

str_file.close()

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


################################################################################
# test                                                                         #
################################################################################

net=Snn()
net.read_net(net_file_name)

avg_error=0

avg_out     =[0,0,0,0,0,0,0]
state_sums  =[0,0,0,0,0,0,0]
avg_real_out=[0,0,0,0,0,0,0]

for dom_id in dom_ids:
    for i in train_index[dom_id]:
        net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
        avg_error += net.ce_error(states[dom_id][i])

        state=states_int[dom_id][i]
        state_sums[state] += 1

        # for k in range(7):
        #     print(net.act_o[k],end=' ')
        # print(state)

        avg_real_out[state] += net.act_o[state]
        for j in range(7):
            avg_out[j] += net.act_o[j]

avg_error /= train_res_sum
for i in range(7):
    avg_out[i] /= train_res_sum
    avg_real_out[i] /=state_sums[i]

print('avg_err',avg_error,'\n')
print('avg_out')
for i in range(7):
    print(avg_out[i])

print('\navg_real_out')
for i in range(7):
    print(avg_real_out[i])

################################################################################
# end                                                                          #
################################################################################
