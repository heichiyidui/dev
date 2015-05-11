#!/usr/bin/env python3
dom_ids=open('index/s.ls').read().split()

s_file = open('../cath/index/cath_s35.cseq')
dom_lens={}
for line in s_file:
    s_line = s_file.readline()
    dom_id = line.strip()[1:]
    if dom_id not in dom_ids :continue 
    
    c_seq = s_line.strip()
    dom_len=0
    for i in range(len(c_seq)):
        if c_seq[i].isupper():
            dom_len += 1
    
    dom_lens[dom_id]=dom_len 
s_file.close()

c_file=open('../cath/index/cath_s35.condef')
for line in c_file:
    c_line=c_file.readline().strip()
    dom_id=line.strip()[1:]
    if dom_id not in dom_ids: continue 
    
    global_contacts=[]
    for col in c_line.split():
        c_i=int(col.split('-')[0])
        c_j=int(col.split('-')[1])
        if c_j-c_i <= 4 : continue 
        global_contacts.append([c_i,c_j])
    
    print('>'+dom_id,dom_lens[dom_id],len(global_contacts))
    for cont in global_contacts:
        print(cont[0],cont[1],end=' ')
    print()
    
c_file.close()
