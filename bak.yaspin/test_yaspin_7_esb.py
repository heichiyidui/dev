#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

################################################################################
# read the data set                                                            #
################################################################################

import sys
dom_ids = open(sys.argv[1]).read().split()
# dom_ids = open('index/test.ls').read().split()
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
# read secondary structure 

import dssp_to_seven

dom_states={}
dom_targets={}
str_file=open('../cath/index/cath_s35.dssp') 
for line in str_file:
    dom_id = line.strip()[1:]
    ss_line = str_file.readline().strip()
    if dom_id not in dom_ids: continue
    
    ss_states = dssp_to_seven.dssp_to_seven(ss_line)
    states=[]
    targets=[]
    for ss in ss_states:
        if ss == -1:
            states.append('X')
            targets.append('X')
            continue
        states.append(ss)
        target=[0.0]*7
        target[ss]=1.0
        targets.append(target)

    dom_targets[dom_id]=targets 
    dom_states[dom_id]=states

str_file.close()

################################################################################
# get the emission probability profile                                         #
################################################################################

import glob
from snn import Snn
net_file_names=glob.glob('res/y7_nets/y_7*.net')
nets=[]
for file_name in net_file_names:
    net=Snn()
    net.read_net(file_name)
    nets.append(net)

from numpy import *
back_ground=array([0.41084,0.03981,0.14190,0.03981,0.03882,0.29002,0.03879])

for dom_id in dom_ids:
    emmision_p=[]
    for i in range(0,len(pssms[dom_id])-IN_SIZE+20,20):
        out = array([0.0] *7)
        for net in nets:
            net.propagate(pssms[dom_id][i:i+IN_SIZE])
            out+=net.act_o

        out/=back_ground
        out/=sum(out)
        emmision_p.append(out)

################################################################################
# forward pass                                                                 #
################################################################################
    for_trans_p=array([\
[0.78929107,0.09731962,0.02223813,0.00000004,0.09107164,0.00007945,0.00000004],\
[0.00000422,0.00000041,0.86055648,0.13943510,0.00000040,0.00000298,0.00000040],\
[0.06049940,0.00000012,0.69271099,0.24179074,0.00499493,0.00000373,0.00000011],\
[0.95370564,0.00000041,0.00000146,0.00000041,0.04625785,0.00003383,0.00000040],\
[0.00000433,0.00000042,0.00000150,0.00000042,0.00000041,0.99999251,0.00000041],\
[0.00012207,0.00000712,0.00000161,0.00000006,0.00000005,0.86599728,0.13387180],\
[0.96188940,0.02721482,0.01089148,0.00000042,0.00000041,0.00000306,0.00000041]])
    t_for_trans_p=for_trans_p.transpose()

    forward_p=[]
    f0_p = for_trans_p[0]*emmision_p[0]
    f0_p/= sum(f0_p)
    forward_p.append(f0_p)

    for i in range(1,len(emmision_p)):
        fi_p=[]
        for j in range(7):
            p = sum(t_for_trans_p[j]*forward_p[i-1]) * emmision_p[i][j]
            fi_p.append(p)
        fi_p /= sum(fi_p)
        forward_p.append(fi_p)

################################################################################
# backward pass                                                                #
################################################################################
    back_emmision_p=[]
    for i in range(len(emmision_p)):
        back_emmision_p.append(emmision_p[-1-i])

    back_trans_p=array([\
[0.78932110,0.00000004,0.02152862,0.09535469,0.00000004,0.00008871,0.09370680],\
[0.97342516,0.00000041,0.00000146,0.00000041,0.00000040,0.00005439,0.02651777],\
[0.06249721,0.24178836,0.69272900,0.00000012,0.00000011,0.00000373,0.00298148],\
[0.00000422,0.13943080,0.86056078,0.00000041,0.00000040,0.00000298,0.00000040],\
[0.93431352,0.00000042,0.01823547,0.04744671,0.00000041,0.00000306,0.00000041],\
[0.00010935,0.00000006,0.00000161,0.00000429,0.13394362,0.86594101,0.00000005],\
[0.00000434,0.00000042,0.00000150,0.00000042,0.00000041,0.99999251,0.00000041]])
    t_back_trans_p=back_trans_p.transpose()

    backward_p=[]
    b0_p = back_trans_p[0]*back_emmision_p[0]
    b0_p/= sum(b0_p)
    backward_p.append(b0_p)

    for i in range(1,len(back_emmision_p)):
        bi_p=[]
        for j in range(7):
            p = sum(t_back_trans_p[j]*backward_p[i-1]) * back_emmision_p[i][j]
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
