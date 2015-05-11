#!/usr/bin/env python3

aa_code_to_int = {'A': 0 ,'R': 1 ,'N': 2 ,'D': 3 ,'C': 4 ,\
                  'Q': 5 ,'E': 6 ,'G': 7 ,'H': 8 ,'I': 9 ,\
                  'L': 10,'K': 11,'M': 12,'F': 13,'P': 14,\
                  'S': 15,'T': 16,'W': 17,'Y': 18,'V': 19,};

aln_file_names=open('t.ls').read().split()

import numpy
for aln_file_name in aln_file_names:
    aln=open('seq_aln/'+aln_file_name).read().split()
    
    dom_seq=aln[0]
    dom_len=len(aln[0])
    
    ofile=open('n_mat/'+aln_file_name,'w')
    for seq in aln[1:]:
        n_matrix=numpy.zeros((20,20),int )
        for i in range(dom_len):
            if seq[i] == '-' : 
                continue
            a1=aa_code_to_int[seq[i]]
            a2=aa_code_to_int[dom_seq[i]]
            n_matrix[a1][a2] += 1
            n_matrix[a2][a1] += 1
        for i in range(20):
            for j in range(20):
                ofile.write(str(n_matrix[i][j]) + ' ')
        ofile.write('\n')
    ofile.close()
    
