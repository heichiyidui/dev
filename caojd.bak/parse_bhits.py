#!/usr/bin/env python3
import sys
ifile_names=[sys.argv[1]]

for ifile_name in ifile_names:
    ifile=open(ifile_name)
    
    dom_id=''
    gi_ids=[]
    gi_start={}
    gi_end={}
    gi_seq_id={}
    # use one segment for each gi id 
    
    for line in ifile:
        if line.startswith('# BLAST'): # done, output 
            if dom_id != '':
                for gi_id in gi_ids:
                    print(dom_id,gi_id,\
                          gi_start[gi_id],gi_end[gi_id],gi_seq_id[gi_id])
                dom_id=''
                gi_ids=[]
                gi_start={}
                gi_end={}
                gi_seq_id={}
        elif line.startswith('# Query:'):
            dom_id = line.split()[2]
            ifile.readline()
            ifile.readline()
            line=ifile.readline()
            num_hits=int(line.split()[1])
            for i in range(num_hits):
                line=ifile.readline()
                cols=line.split()
                
                gi_id=cols[0].split('|')[1]
                t_gi_start=int(cols[1])
                t_gi_end  =int(cols[2])
                t_seq_id  =float(cols[4])
                
                if gi_id in gi_ids and gi_seq_id[gi_id] > t_seq_id:
                    continue 
                
                gi_start [gi_id]=t_gi_start
                gi_end   [gi_id]=t_gi_end
                gi_seq_id[gi_id]=t_seq_id
                
                if gi_id not in gi_ids:
                    gi_ids.append(gi_id)
            
    # break; 

