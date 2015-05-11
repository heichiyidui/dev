#!/usr/bin/env python3

################################################################################
# read the data set                                                            #
################################################################################

dom_ids = open('index/s.ls').read().split()
while '' in dom_ids: dom_ids.remove('')

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

for dom_id in dom_ids:
    states=[]
    for i in range(len(ss_states[dom_id])):
        ss  = ss_states[dom_id][i]
        expo= expos[dom_id][i]
        if ss == 'X' or expo == 'NA':
            states.append('X')
            continue
        state=ss
        if expo < 0.5:
            state +=3
        states.append(state)
    dom_states [dom_id] = states
    
################################################################################
# get the third order transition table                                         #
################################################################################

from numpy import * 

back_ground=array([0.27088, 0.07150, 0.18948, 0.13997, 0.15003, 0.17815])

trans_table=[]
#######################################
# for the first-order HMM
for i in range(6):
    trans_table.append(([0.0]*6)+back_ground)

for dom_id in dom_ids:
    states=dom_states[dom_id]
    # states.reverse()
    for i in range(len(states)-1):
        if 'X' in states[i:i+2]: continue
        state_in  = states[i] 
        state_out = states[i+1]
        trans_table[state_in][state_out]+=1

for i in range(6):
    trans_table[i]/=sum(trans_table[i])
    print('[',end='')
    for j in range(6):
        print('{:11.9f}'.format(trans_table[i][j]),end=',')
    print('],\\')

################################
# for the third-order model 

# for i in range(216):
#     trans_table.append(([0.0]*6)+back_ground)
#
# for dom_id in dom_ids:
#     states=dom_states[dom_id]
#     # states.reverse()
#     for i in range(len(states)-3):
#         if 'X' in states[i:i+4]: continue
#         state_in  = states[i] + 6*states[i+1] + 36 * states[i+2]
#         state_out = states[i+3]
#         trans_table[state_in][state_out]+=1
#
#         # print(state_in,state_out)
#
# for i in range(216):
#     # trans_table[i]/=sum(trans_table[i])
#     for j in range(6):
#         print(trans_table[i][j],end=' ')
#     print()
#
# #     print('[',end='')
# #     for j in range(6):
# #         print('{:11.9f}'.format(trans_table[i][j]),end=',')
# #     print('],\\')
# # print()

################################################################################
# the end                                                                      #
################################################################################
