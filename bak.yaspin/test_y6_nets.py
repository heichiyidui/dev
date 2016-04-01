#!/usr/bin/env python3

################################################################################
# read the data set                                                            #
################################################################################

dom_ids = open('index/test.ls').read().split()
while '' in dom_ids: dom_ids.remove('')

#######################################
# read PSSM 
WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

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
            expo.append('NA')
        else :
            expo.append(float(exp_str))
    
    expos[dom_id]=expo 

exp_file.close()

#######################################
# read secondary structure 
ss_states  = {}

MAP_8_TO_3={'X':'X',
            'E': 1 ,'B': 1 ,
            'H': 2 ,'G': 2 ,
            'C': 0 ,'S': 0 ,'T': 0 ,'I': 0 }

str_file=open('../cath/index/cath_s35.dssp') 
for line in str_file:
    dom_id = line.strip()[1:]
    ss_line = str_file.readline().strip()
    if dom_id not in dom_ids: continue
    
    ss=[]
    for c in ss_line:
        ss.append(MAP_8_TO_3[c])
    
    ss_states[dom_id]=ss
    
str_file.close()

#######################################
# get residue states
dom_states={}
dom_targets={}
for dom_id in dom_ids:
    states=[]
    targets=[]
    for i in range(len(ss_states[dom_id])):
        ss  = ss_states[dom_id][i]
        expo= expos[dom_id][i]
        if ss == 'X' or expo == 'NA':
            states.append('X')
            targets.append('X')
            continue
        state=ss
        if expo < 0.5:
            state +=3
        states.append(state)
        target=[0.0] * 6
        target[state]=1.0
        targets.append(target)
    dom_states [dom_id] = states
    dom_targets[dom_id] = targets

################################################################################
# test                                                                         #
################################################################################

import sys
from snn import Snn 
from numpy import * 

# net=Snn()
# net.read_net('res/y6_nets/y_6_80_90.net')

import glob
net_file_names=glob.glob('res/y6_nets/*.net')
nets=[]
for file_name in net_file_names:
    net=Snn()
    net.read_net(file_name)
    nets.append(net)

avg_error = 0
avg_out      = array([0.0] * 6)
avg_real_out = array([0.0] * 6)
state_sums = array([0] * 6)

for dom_id in dom_ids:
    for i in range(len(dom_states[dom_id])):
        state=dom_states[dom_id][i]
        if state == 'X' : continue

        state_sums[state] +=1 

        out=array([0.0]*6)
        # net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
        # out=net.act_o
        
        for net in nets:
            net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
            out += net.act_o
        out /= len(nets)

        for j in range(6):
            print(out[j],end=' ')
        print(state)

        avg_error -= dot(dom_targets[dom_id][i],log(out))
        avg_out += out
        avg_real_out[state] += out[state]

# res_sum = sum(state_sums)
# avg_error /= res_sum
# avg_out   /= res_sum
# avg_real_out /= state_sums

# print('avg_error',avg_error)
# print('avg_out')
# for i in range(6):
#     print(avg_out[i])
# print('avg_real_out')
# for i in range(6):
#     print(avg_real_out[i])

################################################################################
# the end                                                                      #
################################################################################
