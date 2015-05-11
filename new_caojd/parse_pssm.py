#!/usr/bin/env python3

dom_ids=open('t.ls').read().split()

dom_seqs={}
ifile=open('index/cath_s35.seq')
for line in ifile:
    dom_id=line[1:-1]
    seq_line=ifile.readline().strip()
    if dom_id not in dom_ids:
        continue
    dom_seqs[dom_id]=seq_line
ifile.close()


for dom_id in dom_ids:
    ifile=open('pssm/'+dom_id)
    
    for i in range(3):
        line=ifile.readline()
    v_loc=line.find('V')
    
    out_line=''
    for line in ifile:
        if line=='\n':
            break;
        res_num=int(line[0:5])-1
        aa = line[6]
        if dom_seqs[dom_id][res_num] != aa:
            print('error seq matching!',dom_id,res_num)
        out_line+=line[9:v_loc+1]+' '
        
    ifile.close()
    print('>'+dom_id)
    print(' '.join(out_line.split()))
    
    if len(out_line.split())/20 != len(dom_seqs[dom_id]):
        print('error pssm matching!',dom_id,res_num)
    