#!/usr/bin/env python3

sap_file_names=open('t.ls').read().split()

# read consensus sequences 
domain_ids=set()
for pair_id in sap_file_names:
    domain_ids.add(pair_id.split('_')[0])
    domain_ids.add(pair_id.split('_')[1])

cseqs={}
cseq_num_maps={}

for dom_id in domain_ids:
    cseq=''.join(open('cseq/'+dom_id).read().split('\n')[1:])
    cseqs[dom_id]=cseq

    dom_num_2_seq_num={}
    j=0
    for i in range(len(cseq)):
        if cseq[i].isupper():
            dom_num_2_seq_num[j]=i
            j+=1
    cseq_num_maps[dom_id]=dom_num_2_seq_num

# read SAP alignments

for sap_file_name in sap_file_names:
    dom_id1=sap_file_name.split('_')[0]
    dom_id2=sap_file_name.split('_')[1]
    
    ifile=open('sap_aln/'+sap_file_name)
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
    
    # the output lines to be added into the multiple alignment files 
    alned_dom1 = ['-'] * len(cseqs[dom_id2])
    alned_dom2 = ['-'] * len(cseqs[dom_id1])

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

        res_num1 = cseq_num_maps[dom_id1][res_num1]
        res_num2 = cseq_num_maps[dom_id2][res_num2]

        if AA2 != cseqs[dom_id2][res_num2] or AA1 != cseqs[dom_id1][res_num1] :
            print(line,'wrong AA or residue number',sap_file_name)

        if float(cols[2]) < 2:
            continue 

        alned_dom1[res_num2]=AA1
        alned_dom2[res_num1]=AA2

    # write to the multiple alignment files 
    if alned_dom2.count('-') < len(alned_dom2):
        ofile=open('seq_aln/'+dom_id1,'a')
        ofile.write(''.join(alned_dom2)+'\n')
        ofile.close()

    if alned_dom1.count('-') < len(alned_dom1):
        ofile=open('seq_aln/'+dom_id2,'a')
        ofile.write(''.join(alned_dom1)+'\n')
        ofile.close()

