#!/usr/bin/env python3

aa_code_to_int = {'A': 0 ,'R': 1 ,'N': 2 ,'D': 3 ,'C': 4 ,\
                  'Q': 5 ,'E': 6 ,'G': 7 ,'H': 8 ,'I': 9 ,\
                  'L': 10,'K': 11,'M': 12,'F': 13,'P': 14,\
                  'S': 15,'T': 16,'W': 17,'Y': 18,'V': 19,};

ids = open('t.ls').read().split()

con_defs=[]
for d_id in ids:
    con_def=[]
    ifile = open('conDef/'+d_id)
    for line in ifile:
        cols=line.split()
        con_def.append( (int(cols[0]),int(cols[1])) )
    ifile.close()
    con_defs.append(con_def)

################################################################################
# to get the Pi 
pi_mat =  [0] * 400 
for i in range(len(ids)):
    aln = open('seq_aln/'+ids[i]).read().split()
    dom_seq=aln[0]
    for (loc_i,loc_j) in con_defs[i]:
        if loc_i + 5 > loc_j :
            continue
        a1 = aa_code_to_int[dom_seq[loc_i]]
        a2 = aa_code_to_int[dom_seq[loc_j]]
        pi_mat[a1*20 + a2] += 1

for c in pi_mat:
    print( c/sum(pi_mat) )

#for i in range(20):
#    for j in range(20):
#        print(pi_mat[i*20+j]/sum(pi_mat),end=' ')
#    print('')

