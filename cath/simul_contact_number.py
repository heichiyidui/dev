#!/usr/bin/env python3

# to simulate how many contacts a sequence should have

#######################################
# contact number counts with residue type
residue_contact_counts = {
    'A':[24,441,3025,10736,24163,38861,42455,33307,19838,8553,2693,680,98,11],
    'C':[2,56,253,786,2176,4156,6086,6895,5241,2931,1265,424,118,17,2],
    'D':[84,1040,5339,15408,25537,29803,23308,14716,8306,4007,1651,577,\
         123,27,6],
    'E':[149,1730,7766,19508,31967,34163,25795,15606,8795,4379,2045,766,\
         215,44,5,1],
    'F':[20,329,1360,3049,5334,7232,9117,11754,13903,14307,11646,7084,\
         3274,991,217,29,2,1],
    'G':[261,3388,13011,26216,38475,38277,25499,12411,4820,1470,377,76,11,1],
    'H':[53,586,2203,4681,7258,8992,8650,7567,5714,3672,1933,774,219,63,5,0,1],
    'I':[10,165,945,2875,6390,10989,15435,19397,23720,22888,15907,7748,\
         2840,789,126,15,0,1],
    'K':[196,1852,7842,18191,27763,28556,21433,12742,6373,2621,922,353,\
         113,24,0,1],
    'L':[20,408,1901,5453,11484,18232,23857,29850,36033,35568,26680,14568,\
         5715,1538,266,30,2],
    'M':[39,379,1000,2241,3667,5036,5734,6602,7228,6449,4443,2212,757,162,34,4],
    'N':[91,1036,4378,10911,16794,19966,16930,11192,6669,3358,1299,428,102,\
         10,3],
    'P':[11,567,4422,12216,19430,20674,17159,12857,8295,4118,1502,397,85,6,0,1],
    'Q':[82,730,3377,8580,15118,17859,14976,9901,6097,3271,1437,493,121,29,3],
    'R':[304,1854,6179,13983,20761,22555,19931,14462,8855,4565,1959,776,234,\
         52,7,1],
    'S':[27,560,3987,13166,25041,30828,26838,17074,9462,4090,1432,394,77,10,],
    'T':[21,360,2440,9254,19117,25376,23257,17105,11916,6954,3208,1092,266,\
         60,7],
    'V':[9,263,1477,4700,10998,18672,24925,29638,30713,22801,12069,4642,1465,\
         307,50,6],
    'W':[11,141,482,1187,1974,2698,3417,4084,4689,4492,3663,2208,1010,352,83,\
         8,1],
    'Y':[38,408,1604,3722,6502,8633,10718,12272,11599,9331,6296,3446,1389,\
         384,94,11,2,1]
                         }

# for random choosing of a contact number, given the residue type
residue_contact_number_distribution={}
for AA in residue_contact_counts.keys():
    residue_contact_number_distribution[AA] = [0]*residue_contact_counts[AA][0]

    for i in range(1,len(residue_contact_counts[AA])):
        residue_contact_number_distribution[AA]\
            .extend([i]*residue_contact_counts[AA][i])

# read the sequences
seqs = {}
ifile=open('index/cath_s35.seq')
for line in ifile:
    dom_id = line[1:-1]
    seq = ifile.readline()[:-1]
    seqs[dom_id] = seq
ifile.close()

# read the number of contacts
contact_nums = {}
ifile=open('index/cath_s35.condef')
for line in ifile:
    dom_id = line[1:-1]
    cols = ifile.readline()[:-1].split()
    contact_nums[dom_id]=len(cols)
ifile.close()

#######################################
# simulate possible contact number give the sequence

import numpy as np

for dom_id in seqs.keys():
    simu_cont_sums = []
    for l in range(10):
        simu_cont_sum = 0
        for AA in seqs[dom_id]:
            simu_cont_sum += \
                np.random.choice(residue_contact_number_distribution[AA])
        simu_cont_sum *= 0.5
        simu_cont_sums.append(simu_cont_sum)

    print(dom_id,np.mean(simu_cont_sums),contact_nums[dom_id])
    break;
