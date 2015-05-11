#!/usr/bin/env python3

dom_ids=open('t.ls').read().split()

seqs={}
ifile=open('index/cath_s35.seq')
for line in ifile:
    dom_id=line.strip()[1:]
    seq=ifile.readline().strip()
    if dom_id in dom_ids: seqs[dom_id]=seq 
ifile.close()

stride_ss={}
ifile=open('index/cath_s35.stride')
for line in ifile:
    dom_id=line.strip()[1:]
    ss=ifile.readline().strip()
    if dom_id in dom_ids: stride_ss[dom_id]=ss 
ifile.close()

stride_acc={}
ifile=open('index/cath_s35.stride.acc')
for line in ifile:
    dom_id=line.strip()[1:]
    accs=ifile.readline().strip().split()
    if dom_id in dom_ids: stride_acc[dom_id]=accs
ifile.close()

#######################################
# read the DSSP file 

for dom_id in dom_ids:
    seq=seqs[dom_id]
    acc=stride_acc[dom_id]
    ss=list(stride_ss[dom_id])
    
    ifile=open('dssp_ss/'+dom_id)
    for line in ifile:
        if line.startswith('  #  RESIDUE') :
            break;
    
    for line in ifile:
        if line[13]=='!':
            continue
        seq_index=int(line[5:10])-1
        aa=line[13]
        if aa != seq[seq_index]:
            print('error readling ',dom_id, seq_index)
        
        res_ss=line[16]
        if res_ss==' ' : res_ss='C'
        ss[seq_index]=res_ss
                
        res_acc=line[34:38].strip()
        acc[seq_index]=res_acc
    
    ifile.close()
        
    print('>'+dom_id)
    # print(''.join(ss))
    # print(' '.join(acc))
    
    
    

