#!/usr/bin/env python3

seqs = {}
ifile=open('index/cath_s35.seq')
for line in ifile:
    dom_id = line[1:-1]
    seq = ifile.readline()[:-1]
    seqs[dom_id] = seq
ifile.close()

for dom_id in seqs.keys():
    pssm=[]

    ifile=open('pssm/'+dom_id)
    for i in range(3):
        ifile.readline()

    sum_pssm_aa = 0
    for line in ifile:
        if line=='\n':
            break;
        cols = line.split()

        res_num = int(cols[0]) - 1
        aa = cols[1]

        if seqs[dom_id][res_num] != aa:
            print('mismatch at:',dom_id,res_num)

        pssm_aa = cols[2:22]
        sum_pssm_aa += 1

        pssm.append(' '.join(pssm_aa))
    ifile.close()

    if sum_pssm_aa != len(seqs[dom_id]) :
        print('missing PSSM assignment at:',dom_id)

    print('>'+dom_id)
    print('|'.join(pssm))
