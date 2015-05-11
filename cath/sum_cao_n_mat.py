#!/usr/bin/env python3

aa_code_to_int = {'A': 0 ,'R': 1 ,'N': 2 ,'D': 3 ,'C': 4 ,\
                  'Q': 5 ,'E': 6 ,'G': 7 ,'H': 8 ,'I': 9 ,\
                  'L': 10,'K': 11,'M': 12,'F': 13,'P': 14,\
                  'S': 15,'T': 16,'W': 17,'Y': 18,'V': 19,};

d_ids = open('t.ls').read().split()

################################################################################
# to put the CAO distances into 99 bins/chunks
dises=[]
for d_id in d_ids:
    cao_dises = open('cao_dis/'+d_id).read().split() 
    
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
    sum_n_mats.append([0]*160000)

# n_mat[0] is the diagonal matrix of Pi 
Pi = open('Pi').read().split()
for i in range(400):
    for j in range(400):
        if i!=j :
            continue
        sum_n_mats[0][i*400+j] = float(Pi[i])

for d_id in d_ids:
    cao_dises = open('cao_dis/'+d_id).read().split() 
    
    # read the contact definitions 
    con_def=[]
    ifile = open('conDef/'+d_id)
    for line in ifile:
        cols=line.split()
        # using global contacts only 
        if int(cols[0]) + 5 > int(cols[1]):
            continue 
        con_def.append( ( int(cols[0]),int(cols[1]) ) )
    ifile.close()
    
    for i in range(len(cao_dises)):
        cao_dises[i] = float(cao_dises[i])
    
    aln = open('seq_aln/'+d_id).read().split()
    dom_seq = aln[0]
    
    for i in range(len(cao_dises)):
        seq = aln[i+1]
        n_mat = [0]*160000
        for (loc_i,loc_j) in con_def:
            if seq[loc_i] == '-' or seq[loc_j] =='-' :
                continue 
            a0_1 = aa_code_to_int[dom_seq[loc_i]]
            a0_2 = aa_code_to_int[dom_seq[loc_j]]
            c0 = a0_1 * 20 + a0_2
            
            a1_1 = aa_code_to_int[seq[loc_i]]
            a1_2 = aa_code_to_int[seq[loc_j]]
            c1 = a1_1 * 20 + a1_2 
            
            n_mat[c0 * 400 + c1] += 1
            n_mat[c1 * 400 + c0] += 1
        
        chunk_id = dis_to_chunk[cao_dises[i]]
        for l in range(160000):
            sum_n_mats[chunk_id][l] += n_mat[l]
    
#    print(d_id)
#    break;

for i in range(100):
    for j in range(160000):
        print(sum_n_mats[i][j],end=' ')
    print('')

