#!/usr/bin/env python3

domain_seq={}
ifile=open('index/cath_s35.seq')
for line in ifile:
    dom_id = line[1:-1]
    domain_seq[dom_id]=ifile.readline()[:-1]
ifile.close()

import  os
for sap_file_name in os.listdir("./sap_aln"):
    dom_id1=sap_file_name.split('_')[0]
    dom_id2=sap_file_name.split('_')[1]
    ifile=open('./sap_aln/'+sap_file_name)

    for line in ifile:
        if line.startswith(' score = '):
            sap_score = float(line.split()[2])
            break;

    if sap_score < 200:
        ifile.close()
        continue

    ifile.readline()
    ifile.readline() # The Percent sel on aln line
    ifile.readline() # The Percent aln in sel line

    out_cols = [] # the output for this alignment

    for line in ifile:
        if line.startswith('dompdb/'):
            break;

        cols=line.split()

        if len(cols) != 5:
            print(line,'line misformatted',sap_file_name)
            break;

        AA2      = cols[0][-1]
        res_num2 = int(cols[1]) - 1
        AA1      = cols[4][0]
        res_num1 = int(cols[3]) - 1
        if AA2 != domain_seq[dom_id2][res_num2]:
            print('WTF!', line, sap_file_name)
        if AA1 != domain_seq[dom_id1][res_num1]:
            print('WTF!', line, sap_file_name)

        if float(cols[2]) < 2.0:
            continue
        out_cols.append('{:d}-{:d}'.format(res_num1,res_num2))
    ifile.close()

    print('>'+sap_file_name)
    print(' '.join(out_cols))

