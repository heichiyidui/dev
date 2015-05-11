#!/usr/bin/env python3

################################################################################
# test networks on predicting DSSP 3->10 states definitions                    #
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
        target=[0.0]*10
        if ss == 'X' or expo == 'X':
            states.append('X')
            targets.append(target)
            continue
        state=MAP_5_TO_i[ss]
        if expo < 0.5:
            state +=5
        states.append(state)
        target[state]=1
        targets.append(target)
    
    dom_states[dom_id]=states 
    dom_targets[dom_id]=targets
    
################################################################################
# get the emission probability profile                                         #
################################################################################
net=Snn()
net.read_net('res/y10_nets/y_10_80_02.net')

from numpy import *
back_ground=array([0.27088,0.07150,0.01666,0.15015,0.02267,\
                   0.13997,0.15003,0.01241,0.15256,0.01318])
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
[[0.578968,0.045590,0.044246,0.000038,0.012187, \
  0.236287,0.052066,0.024058,0.000017,0.006543],\
 [0.126221,0.316984,0.001250,0.000006,0.001800, \
  0.108681,0.435115,0.005831,0.000000,0.004111],\
 [0.000000,0.000000,0.000000,0.802679,0.000000, \
  0.000000,0.000000,0.000000,0.197321,0.000000],\
 [0.083508,0.000076,0.000000,0.453707,0.000628, \
  0.040624,0.000368,0.000000,0.420245,0.000843],\
 [0.094713,0.001592,0.000000,0.004233,0.505418, \
  0.132265,0.015231,0.000000,0.015719,0.230830],\
 [0.465857,0.054522,0.035308,0.000030,0.015338, \
  0.274508,0.106307,0.036053,0.000012,0.012065],\
 [0.091374,0.194696,0.001212,0.000003,0.001676, \
  0.107199,0.595905,0.005026,0.000000,0.002908],\
 [0.000000,0.000000,0.000000,0.595902,0.000000, \
  0.000000,0.000000,0.000000,0.404098,0.000000],\
 [0.040373,0.000247,0.000000,0.398748,0.002612, \
  0.028008,0.001409,0.000000,0.526589,0.002013],\
 [0.123901,0.007520,0.000000,0.021474,0.382579, \
  0.116163,0.051742,0.000000,0.037478,0.259144]])

    t_for_trans_p=for_trans_p.transpose()

    forward_p=[]
    f0_p=for_trans_p[0]*emission_p[0]
    f0_p/=sum(f0_p)
    forward_p.append(f0_p)

    for i in range(1,len(emission_p)):
        fi_p=[]
        for j in range(10):
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
[[0.579076,0.034643,0.000000,0.048167,0.008240, \
  0.247261,0.052671,0.000000,0.023668,0.006275],\
 [0.166130,0.316973,0.000000,0.000161,0.000505, \
  0.105432,0.408885,0.000000,0.000528,0.001388],\
 [0.691079,0.005358,0.000000,0.000000,0.000000, \
  0.292652,0.010912,0.000000,0.000000,0.000000],\
 [0.000065,0.000003,0.089115,0.453708,0.000638, \
  0.000027,0.000003,0.049277,0.405277,0.001885],\
 [0.140052,0.005678,0.000000,0.004159,0.505253, \
  0.093543,0.011103,0.000000,0.017595,0.222618],\
 [0.445158,0.056187,0.000000,0.044137,0.021676, \
  0.274442,0.116393,0.000000,0.030927,0.011081],\
 [0.090346,0.207187,0.000000,0.000369,0.002299, \
  0.097890,0.595930,0.000000,0.001433,0.004546],\
 [0.504504,0.033557,0.000000,0.000000,0.000000, \
  0.401194,0.060745,0.000000,0.000000,0.000000],\
 [0.000030,0.000000,0.021552,0.413431,0.002333, \
  0.000011,0.000000,0.032874,0.526533,0.003237],\
 [0.129259,0.022289,0.000000,0.009606,0.396699, \
  0.126492,0.033107,0.000000,0.023315,0.259233]])

    t_back_trans_p=back_trans_p.transpose()

    backward_p=[]
    b0_p = back_trans_p[0]*back_emission_p[0]
    b0_p/= sum(b0_p)
    backward_p.append(b0_p)

    for i in range(1,len(back_emission_p)):
        bi_p=[]
        for j in range(10):
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
              dom_states[dom_id][i])

################################################################################
# end                                                                          #
################################################################################
