#!/usr/bin/env python3 

WING_SIZE = 8 # 17 residues window as PSSM input 
IN_SIZE = (WING_SIZE * 2 + 1) * 20 
WIN_LENGTH=9 # to add 9 extra AAs to the termini in extracting subject sequences

temp_file_name_base='t_93482168'
Y6NETS='res/y6_nets/y_6*.net'
################################################################################
cons_seq = ''
ifile=open(temp_file_name_base+'.pssm')

pssm  = [-2.0] * 20 * WING_SIZE
t_cons_seq=''
ifile.readline()
line=ifile.readline()
if not line.startswith('Last position-specific scoring matrix'):
    print('This does not look like a PSIBlast PSSM file. Exit',file=sys.stderr)
    exit()
ifile.readline()
for line in ifile:
    if line=='\n':break;
    aa = line[6]
    t_cons_seq+=aa
    for i in range(9,69,3):
        pssm.append(float(line[i:i+3]))
ifile.close()
pssm += [-2.0] * 20 * WING_SIZE

if cons_seq != '':
    if t_cons_seq != cons_seq.upper():
        print('Wrong consensus sequence? Exit',file=sys.stderr)
        exit()
else:
    cons_seq = t_cons_seq

import glob
from snn import Snn
from numpy import *
import sys 
################################################################################
# 7. YASPIN 6-state residue exposure prediction                                #
################################################################################

#######################################
# 7.1  PSSM and consensus sequence are loaded 

#######################################
# 7.2 load the neural networks

print('Loading network files ',Y6NETS,file=sys.stderr)

net_file_names=glob.glob(Y6NETS)
nets=[]
for file_name in net_file_names:
    net=Snn()
    net.read_net(file_name)
    nets.append(net)

#######################################
# 7.3 get the emmision probabilities

back_ground=array([0.27088, 0.07150, 0.18948, 0.13997, 0.15003, 0.17815])

emmision_p=[]
for i in range(0,len(pssm)-IN_SIZE+20,20):
    out = array([0.0] *6)
    for net in nets:
        net.propagate(pssm[i:i+IN_SIZE])
        out+=net.act_o

    out/=back_ground
    out/=sum(out)
    emmision_p.append(out)

#######################################
# 7.4 forward pass 
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

#######################################
# 7.5 backward pass 
back_emmision_p=[]
for i in range(len(emmision_p)):
    back_emmision_p.append(emmision_p[-1-i])

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

#######################################
# 7.6 get marginal likelihood and print the results
marginal_like=[]
for i in range(len(forward_p)):
    m_p=forward_p[i]*backward_p[-i-1]
    m_p/=sum(m_p)
    marginal_like.append(m_p)

y6_lables={0:'E',1:'E',2:'E',3:'B',4:'B',5:'B'}
y6_confi_thres=[0.51783,0.60193,0.68635,0.76675,\
                0.83780,0.89489,0.93645,0.96403,0.98211]
y6_results=[]
for i in range(len(marginal_like)):
    aa_confidence=0
    for thres in y6_confi_thres:
        if thres < marginal_like[i].max():
            aa_confidence += 1
    aa_lable=y6_lables[marginal_like[i].argmax()]
    y6_results.append(aa_lable+' '+str(aa_confidence))

################################################################################

seq_i=0
for i in range(len(y6_results)):
    if not cons_seq[i].isupper():
        continue
    seq_i +=1
    print('{:5d}'.format(seq_i),cons_seq[i],y6_results[i])

    