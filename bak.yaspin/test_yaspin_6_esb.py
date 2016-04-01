#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

################################################################################
# test networks on predicting DSSP 3->6 states definitions                     #
################################################################################

WING_SIZE = 8 # 17 residues windows PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 

################################################################################
# read the data set                                                            #
################################################################################
import sys
dom_ids = open(sys.argv[1]).read().split()
#dom_ids = open('index/test.ls').read().split()
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
# ls res/y10_nets/* > t.ls

import glob
from snn import Snn
net_file_names=glob.glob('res/y6_nets/y_6*.net')
nets=[]
for file_name in net_file_names:
    net=Snn()
    net.read_net(file_name)
    nets.append(net)


from numpy import *
back_ground=array([ 0.27088, 0.07150, 0.18948, 0.13997, 0.15003, 0.17815])

for dom_id in dom_ids:
    emmision_p=[]
    for i in range(0,len(pssms[dom_id])-IN_SIZE+20,20):
        out = array([0.0] *6)
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
    [0.57896799, 0.04558986, 0.05647044, 0.23628717, 0.05206629, 0.03061825],\
    [0.12622133, 0.31698394, 0.00305619, 0.10868119, 0.43511468, 0.00994266],\
    [0.07750024, 0.00025082, 0.49159524, 0.04800465, 0.00211254, 0.38053650],\
    [0.46585700, 0.05452206, 0.05067634, 0.27450794, 0.10630734, 0.04812933],\
    [0.09137447, 0.19469631, 0.00289129, 0.10719873, 0.59590524, 0.00793397],\
    [0.04373904, 0.00076784, 0.41511395, 0.03257784, 0.00503461, 0.50276674]])
    
    t_for_trans_p=for_trans_p.transpose()

    forward_p=[]
    f0_p = for_trans_p[0]*emmision_p[0]
    f0_p/= sum(f0_p)
    forward_p.append(f0_p)

    for i in range(1,len(emmision_p)):
        fi_p=[]
        for j in range(6):
            p = sum(t_for_trans_p[j]*forward_p[i-1]) * emmision_p[i][j]
            fi_p.append(p)
        fi_p /= sum(fi_p)
        forward_p.append(fi_p)

################################################################################
# backward pass                                                                #
################################################################################
    back_emmision_p=[]
    for i in range(len(emmision_p)):
        back_emmision_p.append(emmision_p[-(i+1)])

    back_trans_p=array([\
    [0.57907551, 0.03464334, 0.05640696, 0.24726086, 0.05267084, 0.02994249],\
    [0.16612961, 0.31697304, 0.00066511, 0.10543215, 0.40888502, 0.00191507],\
    [0.07759901, 0.00115245, 0.49157717, 0.03695407, 0.00228976, 0.39042753],\
    [0.44515754, 0.05618680, 0.06581214, 0.27444203, 0.11639330, 0.04200819],\
    [0.09034615, 0.20718728, 0.00266752, 0.09789001, 0.59592965, 0.00597939],\
    [0.04473148, 0.00398604, 0.40455567, 0.03731340, 0.00668018, 0.50273322]])

    t_back_trans_p=back_trans_p.transpose()

    backward_p=[]
    b0_p = back_trans_p[0]*back_emmision_p[0]
    b0_p/= sum(b0_p)
    backward_p.append(b0_p)

    for i in range(1,len(back_emmision_p)):
        bi_p=[]
        for j in range(6):
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
