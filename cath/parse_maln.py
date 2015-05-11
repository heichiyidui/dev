#!/usr/bin/env python3

domain_ids=open('t.ls').read().split()

for dom_id in domain_ids:
    cons_seq=open('cseq/'+dom_id).read().split('\n')
    cseq_title=cons_seq[0][1:]
    cseq=''.join(cons_seq[1:])
    
    ifile=open('c_align/'+dom_id)
    alignment = ifile.read().split('\n>')
    ifile.close()

    alned_cseq=''
    alned_sub_seqs=[]
    for seq in alignment:
        if seq == '':    continue
        if seq[0] =='>': seq=seq[1:]

        # only the first entry of consensus sequence to be used as aln_cseq
        if seq.split('\n')[0]==cseq_title:
            if alned_cseq == '':
                alned_cseq=''.join(seq.split('\n')[1:])
                continue

        alned_sub_seqs.append(''.join(seq.split('\n')[1:]))

    if alned_cseq.replace('-','') != cseq.upper() :
        print('failed to match cseq in maln',dom_id)
        break;

    for seq in alned_sub_seqs:
        if len(seq) != len(alned_cseq):
            print('length error',seq,dom_id)
    
    # remove consensus sequence gap columns in the MSA
    for i in range(len(alned_sub_seqs)):
        seq = ''
        for j in range(len(alned_cseq)):
            if alned_cseq[j] != '-':
                if alned_sub_seqs[i][j] in '-ARNDCQEGHILKMFPSTWYV':
                    seq += alned_sub_seqs[i][j]
                else:
                    seq += '-'
        alned_sub_seqs[i]=seq 
    alned_cseq=cseq 

    # write output file 
    ofile =open('seq_aln/'+dom_id,'w')
    ofile.write(alned_cseq+'\n')
    for seq in alned_sub_seqs:
        ofile.write(seq+'\n')
    ofile.close()

