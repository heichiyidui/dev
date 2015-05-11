#!/usr/bin/env python3 
##!/home/klinbrc/bin/python3
import os
import sys
import subprocess

################################################################################
# global parameters                                                            #
################################################################################

ALIGNER='~/bin/muscle -quiet -in '

################################################################################
# 1. make the output directions                                                #
################################################################################

# the multiple alignments of the sequences
if not os.path.exists('./c_align'):
    subprocess.call('mkdir ./c_align',shell=True)

# the consensus sequence to be used for the next blast 
if not os.path.exists('./cseq'):
    subprocess.call('mkdir ./cseq',shell=True)

################################################################################
# 2. read the list of query sequence ids                                       #
################################################################################

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])

ifile.close()

################################################################################
# 3 multiple alignment of the subject sequences and the query sequence         #
################################################################################

for id in domain_ids:
    
    ###################################
    # 3.1 muscle alignment 
#    if not os.path.isfile('./c_align_in/'+id):
#        sys.stderr.write('Failed to open blast seqs file c_align_in/'+id+'\n')
#        exit(1)
#    
#    print(id,'aligning sequences...')
#    subprocess.call(ALIGNER+'c_align_in/'+id+' > c_align/'+id,shell=True)
    
    ###################################
    # 3.2 get the query seq title
    
    ifile=open('./c_align_in/'+id)
    for line in ifile: 
        if line[0]==">":
            query_seq_header=line[1:-1]
            break;
    
    ifile.close() # the query sequence should always be the first in the input
    
    ###################################
    # 3.3 get the aligned sequences
    
    ifile = open('c_align/'+id)
    alignment = ifile.read().split('\n>')
    ifile.close()
    
    aln_subj_seqs = []
    aln_query_seq = ''
    
    for seq in alignment:
        if seq == '' : continue 
        
        if seq[0]=='>':
            seq = seq [1:]
        
        if seq.split('\n')[0] == query_seq_header:
            aln_query_seq = ''.join(seq.split('\n')[1:]).upper()
        else:
            aln_subj_seqs.append(''.join(seq.split('\n')[1:]).upper())
    
    # Better make sure the title of the query seq is never in the BLAST hits 
    # But when there is no BLAST hits, aln_subj_seqs is empty, cseq will be 
    # identical to the domain sequence. 
    
    ###################################
    # 3.4 replace the terminal gaps of the aligned subject seqs with 'u' 
    
    u_filled_aln_subj_seqs=[]
    for seq in aln_subj_seqs:
        left_u  = 'u' * (len(seq) - len(seq.lstrip('-')))
        right_u = 'u' * (len(seq) - len(seq.rstrip('-')))
        u_filled_seq = ''.join((left_u,seq.strip('-'),right_u))
        u_filled_aln_subj_seqs.append(u_filled_seq)
    
    aln_subj_seqs = u_filled_aln_subj_seqs
    # quite often the subject sequences can be quite short, 
    # with long gaps on their termini 
    # when we try to get a consensus, those gaps should not count. 
    
    ###################################
    # 3.5 to get the consensus 
    
    if aln_query_seq.find('-') == -1:
        cons_seq = aln_query_seq
    else:
        aln_len  = len(aln_query_seq)
        
        # transpose the aln_subj_seqs matrices 
        trans_aln_subj_seqs = \
            [[row[i] for row in aln_subj_seqs] for i in range(aln_len)]
        
        for i in range(aln_len):
            trans_aln_subj_seqs[i] = ''.join(trans_aln_subj_seqs[i])
            trans_aln_subj_seqs[i] = trans_aln_subj_seqs[i].replace('u','')
        
        cons_seq = ''
        for i in range(aln_len):
            if aln_query_seq[i] != '-':
                cons_seq += aln_query_seq[i]
            else:
                aln_depth = len(trans_aln_subj_seqs[i])
                
                is_cons_c_found = False
                for C in 'ARNDCQEGHILKMFPSTWYV':
                    if     trans_aln_subj_seqs[i].count(C) > aln_depth * 0.5 \
                       and aln_depth > 5 :
                       cons_seq += C.lower()
                       is_cons_c_found = True
                       break;
                
                if     is_cons_c_found == False and aln_depth > 5 \
                   and aln_depth > len(aln_subj_seqs) * 0.7 :
                    
                    if trans_aln_subj_seqs[i].count('-') > aln_depth * 0.5:
                        continue
                    else:
                        cons_seq += 'x'
        
    
    ###################################
    # 3.6 write the consensus seq 
    
    ofile=open('cseq/'+id,'w')
    ofile.write('>cons_'+id+'\n')
    for i in range(0,len(cons_seq),80):
        ofile.write(cons_seq[i:i+80]+'\n')
    
    ofile.close()



