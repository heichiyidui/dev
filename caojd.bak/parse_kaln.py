#!/usr/bin/env python3
##!/home/klinbrc/bin/python3
import sys
import fastaio
from collections import Counter

dom_ids=open(sys.argv[1]).read().split()

for dom_id in dom_ids:
    ifile_name='kaln/'+dom_id 
    seqs=fastaio.read_fasta(ifile_name)
    
    seq_set=set([])
    seq_lines=[]
    for seq in seqs:
        seq_line=seq.split('\n')[1]
        seq_line=seq_line.replace('B','N').replace('Z','Q').replace('X','-')\
                         .replace('O','K').replace('U','C')
        if seq_line not in seq_set:
            seq_set.add(seq_line)
            seq_lines.append(seq_line)
    
    aln_length=len(seq_lines[0])
    u_seq_lines=[] # sequences with termini GAPs replaced with 'u'
    for line in seq_lines:
        line=line.lstrip('-')
        line='u'*(aln_length-len(line))+line
        line=line.rstrip('-')
        line=line+'u'*(aln_length-len(line))
        u_seq_lines.append(line)
    
    # keep columns with domain seq AA 
    dom_aa_pos=[]
    for i in range(aln_length):
        if u_seq_lines[0][i] not in 'u-':
            dom_aa_pos.append(i)
    
    dom_len=len(dom_aa_pos)
    ofile=open('temp/'+dom_id,'w')
    for line in u_seq_lines:
        new_line=''
        for pos in dom_aa_pos:
            new_line+=line[pos]
        
        char_counter=Counter(new_line)
        if (char_counter['u'] + char_counter['-']) > (dom_len * 0.9):
            continue
        
        ofile.write(new_line+'\n')
    ofile.close()
    