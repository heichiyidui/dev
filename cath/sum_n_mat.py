#!/usr/bin/env python3

cao_file_names = open('t.ls').read().split()

################################################################################
# to put the CAO distances into 99 bins/chunks
dises=[]
for cao_file_name in cao_file_names:
    cao_dises = open('cao_dis/'+cao_file_name).read().split() 
    
    for i in range(len(cao_dises)):
        dises.append(float(cao_dises[i]))
    
sorted_dises = sorted(set(dises))

dis_to_chunk = {}
chunk_size = len(sorted_dises) / 99 + 1
for i in range(len(sorted_dises)):
    dis_to_chunk[sorted_dises[i]] = int(i / chunk_size) + 1 

################################################################################
# sum n_mats for each chunk
sum_n_mats=[]
for i in range(100):
    sum_n_mats.append([0]*400)
Pi = open('Pi').read().split()
for i in range(20):
    for j in range(20):
        if i!=j :
            continue
        sum_n_mats[0][i*20+j] = float(Pi[i])

for cao_file_name in cao_file_names:
    cao_dises = open('cao_dis/'+cao_file_name).read().split() 
    
    for i in range(len(cao_dises)):
        cao_dises[i] = float(cao_dises[i])
    
    ifile=open('n_mat/'+cao_file_name)
    for cao_dis in cao_dises:
        chunk_id = dis_to_chunk[cao_dis]
        cols = ifile.readline().split()
        for i in range(400):
            sum_n_mats[chunk_id][i] += int(cols[i])
    ifile.close()
    
#    break;

for i in range(100):
    sum_n=0
    for j in range(400):
        sum_n += sum_n_mats[i][j]
    for j in range(400):
        print(sum_n_mats[i][j]/sum_n,end=' ')
    print('')

