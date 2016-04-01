#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

################################################################################
# test the snn neural network on residue exposure prediction                   #
################################################################################

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

import sys
import copy
import random
from snn import Snn

################################################################################
# read the test set 
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

test_res_sum=0
test_index={}
for dom_id in dom_ids:
    t_index=[]
    for i in range(len(expos[dom_id])):
        if expos[dom_id][i] == 'NA': 
            continue
        t_index.append(i)
    
    random.shuffle(t_index)
    test_index[dom_id]=t_index
    test_res_sum += len(t_index)

################################################################################
# testing 

net=Snn(IN_SIZE,20,2,eta=0.001)
net.read_net(sys.argv[1])

correct_sum=0
for dom_id in dom_ids:
    for i in test_index[dom_id]:
        net.propagate(pssms[dom_id][i*20:i*20+IN_SIZE])
        if expos[dom_id][i] > 0.5:
            if net.act_o[0] > 0.5: correct_sum +=1
        else:
            if net.act_o[0] <=0.5: correct_sum +=1

print(sys.argv[1],correct_sum/test_res_sum)
