#!/usr/bin/env python3
import sys

seqs={}
ifile=open('index/cath_s35.seq')
for line in ifile:
    dom_id = line[1:-1]
    seq = ifile.readline()[:-1]
    seqs[dom_id] = seq
ifile.close()

ALLOWED_AA_SET = set('ARNDCQEGHILKMFPSTWYV')

for dom_id in seqs.keys():
    seq = seqs[dom_id]
    dom_len = len(seq)
    ifile=open('bl_out/'+dom_id)

    for i in range(4):
        ifile.readline()
    cols = ifile.readline().split()
    num_hits = int(cols[1])
    print(num_hits)

    out_subject_seqs = []
    out_subject_ids = []
    for i in range(num_hits):
        cols = ifile.readline()[:-1].split()

        aln_start, aln_end = int(cols[0]) - 1, int(cols[1]) -1
        # It seems aln_end is not needed.

        q_seq, s_seq = cols[2], cols[3]
        subject_id = cols[4]

        if subject_id in out_subject_ids: # found this subject sequence before
            continue

        out_subject_seq = ['-'] * dom_len
        dom_index = aln_start
        for j in range(len(q_seq)):
            if q_seq[j] == '-':
                continue

            if q_seq[j] != seq[dom_index]:
                print('WAT?',dom_id,dom_index,q_seq[j],seq[dom_index])

            out_subject_seq[dom_index] = s_seq[j]
            if s_seq[j] not in ALLOWED_AA_SET:
                out_subject_seq[dom_index] = '-'

            dom_index += 1

        out_subject_seq = ''.join(out_subject_seq)

        if out_subject_seq == seq:
            continue
        if out_subject_seq in out_subject_seqs:
            continue

        diff_aa_sum = 0
        for j in range(dom_len):
            if out_subject_seq[j] == '-':
                continue
            if out_subject_seq[j] == seq[j]:
                continue
            diff_aa_sum += 1
        if diff_aa_sum == 0 :
            continue

        out_subject_ids.append(subject_id)
        out_subject_seqs.append(out_subject_seq)
    ifile.close()

    ofile = open('bl_out2/'+dom_id,'w')
    ofile.write('>'+dom_id+'\n')
    ofile.write(seq+'\n')
    for i in range(len(out_subject_ids)):
        ofile.write('>'+out_subject_ids[i]+'\n')
        ofile.write(out_subject_seqs[i] + '\n')
    ofile.close()
    # break;


