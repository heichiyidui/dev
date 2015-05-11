#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

################################################################################
# test networks on predicting DSSP 3->6 states definitions                     #
################################################################################

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

import sys
import copy
import random
from snn import Snn

################################################################################
# read the data set                                                            #
################################################################################
import sys
dom_ids = open(sys.argv[1]).read().split()
#dom_ids = open('index/t.ls').read().split()
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
# get the emission probability profile                                         #
################################################################################

import glob
net_file_names=glob.glob('res/y6_nets/y_6*.net')
nets=[]
for file_name in net_file_names:
    net=Snn()
    net.read_net(file_name)
    nets.append(net)

from numpy import *
back_ground=array([0.27088, 0.07150, 0.18948, 0.13997, 0.15003, 0.17815])

import third_ord_6_state_trans_table
t_for_trans_p  = array(third_ord_6_state_trans_table.for_trans_p ).transpose()
t_back_trans_p = array(third_ord_6_state_trans_table.back_trans_p).transpose()

for dom_id in dom_ids:
    emission_p=[]
    for i in range(0,len(pssms[dom_id])-IN_SIZE+20,20):
        out = array([0.0] *6)
        for net in nets:
            net.propagate(pssms[dom_id][i:i+IN_SIZE])
            out+=net.act_o

        out/=back_ground
        out/=sum(out)
        emission_p.append(out)
    
################################################################################
# forward pass                                                                 #
################################################################################
    
    forward_p=[]
    forward_p.append(emission_p[0])
    forward_p.append(emission_p[1])
    forward_p.append(emission_p[2])

    for i in range(0,len(emission_p)-3):
        previous_p=array([0.0]*216)
        for l in range(6):
            for m in range(6):
                for n in range(6):
                    int_id = l + 6*m + 36 *n 
                    previous_p[int_id] = forward_p[i  ][l] \
                                        *forward_p[i+1][m] \
                                        *forward_p[i+2][n]
        previous_p /= sum(previous_p)
        fi_p=[]
        for j in range(6):
            p = sum(t_for_trans_p[j]*previous_p) * emission_p[i+3][j]
            fi_p.append(p)
        fi_p /= sum(fi_p)
        forward_p.append(fi_p)

################################################################################
# backward pass                                                                #
################################################################################
    back_emission_p=[]
    for i in range(len(emission_p)):
        back_emission_p.append(emission_p[-1-i])

    backward_p=[]
    backward_p.append(back_emission_p[0])
    backward_p.append(back_emission_p[1])
    backward_p.append(back_emission_p[2])

    for i in range(0,len(back_emission_p)-3):
        previous_p=array([0.0]*216)
        for l in range(6):
            for m in range(6):
                for n in range(6):
                    int_id = l + 6*m + 36 *n 
                    previous_p[int_id] = backward_p[i  ][l] \
                                        *backward_p[i+1][m] \
                                        *backward_p[i+2][n]
        previous_p /= sum(previous_p)

        bi_p=[]
        for j in range(6):
            p = sum(t_back_trans_p[j]*previous_p) * back_emission_p[i+3][j]
            bi_p.append(p)
        bi_p /= sum(bi_p)
        backward_p.append(bi_p)

################################################################################
# summing up forward and backward probabilities                                #
################################################################################
    marginal_like=[]
    for i in range(len(forward_p)):
        m_p = forward_p[i] * backward_p[-i-1]
        m_p/= sum(m_p)
        marginal_like.append(m_p)

    for i in range(len(marginal_like)):
        print(marginal_like[i].argmax(),marginal_like[i].max(),\
              dom_states[dom_id][i])

################################################################################
# end                                                                          #
################################################################################
