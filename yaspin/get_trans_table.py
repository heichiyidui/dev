#!/usr/bin/env python3


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

MAP_5_TO_i={'C':0,'E':1,'Hb':2,'H':3,'G':4}

for dom_id in dom_ids:
    states=[]
    for i in range(len(ss_states[dom_id])):
        ss  = ss_states[dom_id][i]
        expo= expos[dom_id][i]
        if ss == 'X' or expo == 'X':
            states.append('X')
            continue
        state=MAP_5_TO_i[ss]
        if expo < 0.5:
            state +=5
        states.append(state)
    
    dom_states[dom_id]=states 
    

#######################################
# get first order transition table 

# trans_table=[]
# for i in range(10):
#     trans_table.append([0]*10)

# for dom_id in dom_ids:
#     states=dom_states[dom_id]
#     for i in range(len(states)-1):
#         if states[i]  =='X':continue
#         if states[i+1]=='X':continue
#         trans_table[states[i]][states[i+1]] +=1

# for i in range(10):
#     print('[',end='')
#     for j in range(10):
#         print('{:6d}'.format(trans_table[i][j]),end=',')
#     print('],')
# print()

# rev_trans_table=[]
# for i in range(10):
#     rev_trans_table.append([0]*10)

# for dom_id in dom_ids:
#     states=dom_states[dom_id]
#     states.reverse()
#     for i in range(len(states)-1):
#         if states[i]  =='X':continue
#         if states[i+1]=='X':continue
#         rev_trans_table[states[i]][states[i+1]] +=1

# for i in range(10):
#     print('[',end='')
#     for j in range(10):
#         print('{:6d}'.format(rev_trans_table[i][j]),end=',')
#     print('],')
# print()

#######################################
# get third order trans_table table 

trans_table=[]
for i in range(10):
    trans_table.append([0]*1000)
for dom_id in dom_ids:
    states=dom_states[dom_id]
    for i in range(len(states)-3):
        if states[i]  =='X':continue
        if states[i+1]=='X':continue
        if states[i+2]=='X':continue
        if states[i+3]=='X':continue

        s_i = states[i] + states[i+1]*10 + states[i+2] *100
        trans_table[states[i+3]][s_i] +=1

for i in range(10):
    for j in range(1000):
        print('{:6d}'.format(trans_table[i][j]),end=' ')
    print()
print()
print()

rev_trans_table=[]
for i in range(10):
    rev_trans_table.append([0]*1000)
for dom_id in dom_ids:
    states=dom_states[dom_id]
    states.reverse()
    for i in range(len(states)-3):
        if states[i]  =='X':continue
        if states[i+1]=='X':continue
        if states[i+2]=='X':continue
        if states[i+3]=='X':continue

        s_i = states[i] + states[i+1]*10 + states[i+2] *100
        rev_trans_table[states[i+3]][s_i] +=1

for i in range(10):
    for j in range(1000):
        print('{:6d}'.format(rev_trans_table[i][j]),end=' ')
    print()
print()

