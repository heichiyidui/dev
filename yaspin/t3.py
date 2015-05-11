#!/usr/bin/env python3

################################################################################
# read states                                                                  #
################################################################################

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

dom_ids = open('index/t.ls').read().split()
while '' in dom_ids: dom_ids.remove('')

#######################################
# read PSSM 

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
# get states

dom_states={}
dom_targets={}
MAP_5_TO_i={'C':0,'E':1,'Hb':2,'H':3,'G':4}

for dom_id in dom_ids:
    states=[]
    targets=[]
    for i in range(len(ss_states[dom_id])):
        ss  = ss_states[dom_id][i]
        expo= expos[dom_id][i]
        target=[0.0]*2
        if ss == 'X' or expo == 'NA':
            states.append('X')
            targets.append(target)
            continue
        state=0
#         state=MAP_5_TO_i[ss]
        if expo < 0.5:
            state +=1
        states.append(state)
        target[state]=1
        targets.append(target)
    
    dom_states[dom_id]=states 
    dom_targets[dom_id]=targets
    
################################################################################
# test networks                                                                #
################################################################################

import sys
import os
from snn import Snn

net2s=[]
for file_name in os.listdir('res/y2_nets'):
    net=Snn()
    net.read_net('res/y2_nets/'+file_name)
    net2s.append(net)

from numpy import * 

avg_error=0

avg_out     =array([0.0] * 2)
state_sums  =array([0.0] * 2)
avg_real_out=array([0.0] * 2)

test_res_sum=0

for dom_id in dom_ids:
    for i in range(len(dom_states[dom_id])):
        if dom_states[dom_id][i]=='X': continue

        test_res_sum +=1

        out_2=array([0.0,0.0])
        for net in net2s:
            net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
            out_2 += net.act_o
        out_2 /= len(net2s)
        print(out_2,state)

        avg_error -= dot(dom_targets[dom_id][i],log(out_2))

        state=dom_states[dom_id][i]
        state_sums[state] += 1
        print(out_2,state)

        avg_real_out[state] += out_2[state]
        avg_out += out_2

avg_error /= test_res_sum 
avg_out   /= test_res_sum 
avg_real_out /=state_sums

print('avg_err',avg_error,'\n')
print('avg_out',avg_out)
print('avg_real_out',avg_real_out)

################################################################################
# end                                                                          #
################################################################################
