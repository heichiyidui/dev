#!/usr/bin/env python3

################################################################################
# test networks on predicting DSSP 3->7 states definitions                     #
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

dom_ids = open('index/test.ls').read().split()

pssms={}
pssm_file=open('../cath/index/cath_s35.pssm')
for line in pssm_file:
    d_id = line.strip()[1:]
    pssm_line = pssm_file.readline().strip()
    
    if d_id not in dom_ids: continue
    
    pssm  = [-2.0] * 20 * WING_SIZE
    
    for residue_pssm in pssm_line.split():
        for score in residue_pssm.split('|'):
            pssm.append(float(score)) 
    
    pssm += [-2.0] * 20 * WING_SIZE
    pssms[d_id]=pssm 

pssm_file.close()

import dssp_to_seven
states  = {}
states_int={}
str_file=open('../cath/index/cath_s35.dssp') 
for line in str_file:
    d_id = line.strip()[1:]
    ss_line = str_file.readline().strip()
    if d_id not in dom_ids: continue
    
    ss_states = dssp_to_seven.dssp_to_seven(ss_line)
    states_int[d_id]=ss_states
    ss=[]
    for state in ss_states:
        if state == -1: 
            ss.append('X')
            continue
        
        target = [0.0] *7
        target[state]=1.0
        ss.append(target)
    
    states[d_id]=ss 

str_file.close()

################################################################################
# get the emission probability profile                                         #
################################################################################
net=Snn()
net.read_net('res/y7_nets/y_7_80_47.net')

from numpy import *
back_ground=array([0.40952,0.03969,0.14157,0.03969,0.03905,0.29146,0.03902])

for dom_id in dom_ids:
    emission_p=[]
    for i in range(0,len(pssms[dom_id])-IN_SIZE+20,20):
        net.propagate(pssms[dom_id][i:i+IN_SIZE])
        out=net.act_o
        out/=back_ground
        out/=sum(out)
        emission_p.append(out)

################################################################################
# forward pass                                                                 #
################################################################################

    for_trans_p=array(\
   [[0.794441, 0.093144, 0.021889, 0.000008, 0.090494, 0.000016, 0.000008 ],\
    [0.000083, 0.000083, 0.856335, 0.143249, 0.000083, 0.000083, 0.000083 ],\
    [0.058703, 0.000022, 0.711203, 0.225093, 0.004936, 0.000022, 0.000022 ],\
    [0.949564, 0.000083, 0.000083, 0.000083, 0.050021, 0.000083, 0.000083 ],\
    [0.000082, 0.000082, 0.000082, 0.000082, 0.000082, 0.999508, 0.000082 ],\
    [0.000061, 0.000010, 0.000010, 0.000010, 0.000010, 0.875299, 0.124599 ],\
    [0.960154, 0.026318, 0.013200, 0.000082, 0.000082, 0.000082, 0.000082 ]])

    t_for_trans_p=for_trans_p.transpose()

    forward_p=[]
    f0_p=for_trans_p[0]*emission_p[0]
    f0_p/=sum(f0_p)
    forward_p.append(f0_p)

    for i in range(1,len(emission_p)):
        fi_p=[]
        for j in range(7):
            p = sum(t_for_trans_p[j]*forward_p[i-1]) * emission_p[i][j]
            fi_p.append(p)
        fi_p /= sum(fi_p)
        forward_p.append(fi_p)

################################################################################
# backward pass                                                                #
################################################################################
    back_emission_p=[]
    for i in range(len(emission_p)):
        back_emission_p.append(emission_p[-(i+1)])

    back_trans_p=array(\
   [[0.794378, 0.000008, 0.021394, 0.090956, 0.000008, 0.000048, 0.093209 ],\
    [0.972897, 0.000083, 0.000083, 0.000083, 0.000083, 0.000083, 0.026688 ],\
    [0.060070, 0.224989, 0.711359, 0.000022, 0.000022, 0.000022, 0.003517 ],\
    [0.000083, 0.143166, 0.856419, 0.000083, 0.000083, 0.000083, 0.000083 ],\
    [0.931809, 0.000082, 0.018523, 0.049340, 0.000082, 0.000082, 0.000082 ],\
    [0.000020, 0.000010, 0.000010, 0.000010, 0.124640, 0.875299, 0.000010 ],\
    [0.000082, 0.000082, 0.000082, 0.000082, 0.000082, 0.999508, 0.000082 ]])

    t_back_trans_p=back_trans_p.transpose()

    backward_p=[]
    b0_p = back_trans_p[0]*back_emission_p[0]
    b0_p/= sum(b0_p)
    backward_p.append(b0_p)

    for i in range(1,len(back_emission_p)):
        bi_p=[]
        for j in range(7):
            p = sum(t_back_trans_p[j]*backward_p[i-1]) * back_emission_p[i][j]
            bi_p.append(p)
        bi_p /= sum(bi_p)
        backward_p.append(bi_p)

################################################################################
# summing up forward and backward probabilities                                #
################################################################################
    marginal_like=[]
    for i in range(len(forward_p)):
        m_p=forward_p[i]*backward_p[-i-1]
        m_p/=sum(m_p)
        marginal_like.append(m_p)

    for i in range(len(marginal_like)):
        print(marginal_like[i].argmax(),marginal_like[i].max(),\
              states_int[dom_id][i])


################################################################################
# end                                                                          #
################################################################################
