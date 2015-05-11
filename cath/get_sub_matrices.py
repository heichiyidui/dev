#!/usr/bin/env python3

GLOBAL_DIS = 4  # the distance between residues in global contacts is at least 4
PATTEN_WIDTH = 5 # 5 x 5 sub-matrices

########################################
# to get H class domain ids 
domain_ids=[]
CATH_set=set()
ifile=open('index/CathDomainList.S35')
for line in ifile:
    dom_id=line.split()[0]
    CATH='_'.join(line.split()[1:5])
    if CATH not in CATH_set:
        CATH_set.add(CATH)
        domain_ids.append(dom_id)

########################################
# read domain DSSP 
EIGHT_TO_THREE = \
    {'H':'H','G':'H', 'E':'E','B':'E', 'S':'C','T':'C','I':'C','C':'C','X':'C'}
# the last 'X':'C' in the dictionary means the 110 residues with missing 
# DSSP assignment will be treated as 'C'

dom_dssp={}
dom_lens={}
ifile=open('index/cath_s35.dssp')
for line in ifile:
    dom_id=line.strip()[1:]
    dssp=ifile.readline().strip()

    if dom_id not in domain_ids: continue;
    dssp_3=''
    for c in dssp:
        dssp_3 += EIGHT_TO_THREE[c]

    dom_dssp[dom_id] = dssp_3
    dom_lens[dom_id] = len(dssp_3)
ifile.close()

########################################
# read domain contact matrices, write sub-matrices
ss_code = {'C':'C ','H':'H ','E':'E '}
cent_shift = int(PATTEN_WIDTH/2)+1

ifile=open('index/cath_s35.condef')
for line in ifile:
    dom_id = line.strip()[1:]
    cont_line = ifile.readline()

    if dom_id not in domain_ids:
        continue
    
    dom_len = dom_lens[dom_id]

    cont_matrix=[]
    for i in range(dom_len):
        cont_matrix.append([0]*dom_len)

    contacts = cont_line.split()
    for contact in contacts:
        i,j = int(contact.split('-')[0]),int(contact.split('-')[1])
        cont_matrix[i][j] = 1
        cont_matrix[j][i] = 1
    
    for i in range(dom_len - PATTEN_WIDTH):
        for j in range(i,dom_len - PATTEN_WIDTH):

            # global contacts only
            if j-i < GLOBAL_DIS + PATTEN_WIDTH: continue
            
            # remove all zero sub-matrices
            sub_matrix_sum=0
            for m_i in range(PATTEN_WIDTH):
                for m_j in range(PATTEN_WIDTH):
                    sub_matrix_sum += cont_matrix[i+m_i][j+m_j]
            if sub_matrix_sum == 0: continue

            # secondary structure
            i_cent = i + cent_shift
            j_cent = j + cent_shift
            
            ss_i = dom_dssp[dom_id][i_cent]
            ss_j = dom_dssp[dom_id][j_cent]
            
            ss_code_i = ss_code[ss_i]
            ss_code_j = ss_code[ss_j]
            
            print(ss_code_i,ss_code_j,end=' ')

            # print the sub-matrix
            for m_i in range(PATTEN_WIDTH):
                for m_j in range(PATTEN_WIDTH):
                    print(cont_matrix[i+m_i][j+m_j],end=' ')
            print()

ifile.close()

