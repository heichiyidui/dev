#!/usr/bin/env python3 
import sys 
import fastaio

ifile=open(sys.argv[1])
for line in ifile:
    if line.startswith('# BLAST processed'): 
        break # "BLAST processed 5 queries" etc. end of file 
    
    if not line.startswith('# BLASTP 2.2'):
        print('WAT? ')
        print(sys.argv[1],line[:-1])
        break
    
    dom_id = ifile.readline()[:-1].split()[-1]
    dom_seq = fastaio.read_fasta('seq/'+dom_id)[0].split('\n')[1]
    dom_len = len(dom_seq)
        
    ifile.readline() # the RID line 
    ifile.readline() # the nr database line 
    ifile.readline() # the fields line 
    num_hits = int(ifile.readline().split()[1])
    
    out_alns=[]
    for i in range(num_hits):
        cols = ifile.readline()[:-1].split('\t')
        
        aln_start = int(cols[0])
        aln_end   = int(cols[1])
        aln_query = cols[2]
        aln_hit   = cols[3]
        
        pairwise_aln = '.' * (aln_start - 1)
        
        # if len(aln_query) != len(aln_hit) :
        #     print('WTF? Aligned query and subject of different length?')
        #     print(sys.argv[1],cols)
        # if aln_end - aln_start + 1  > dom_len:
        #     print('WTF? aligned part larger than domain length?')
        #     print(sys.argv[1],cols)
        
        aln_index = aln_start - 1
        for i in range(len(aln_query)):
            if aln_query[i] == '-':
                continue 
            
            if aln_query[i] != dom_seq[aln_index]:
                print('WTF? Aligned query changed by BLAST?')
                print(sys.argv[1],cols)
                break;
            
            pairwise_aln += aln_hit[i]
            aln_index += 1
        
        pairwise_aln += '.' * (dom_len - aln_end)
        
        # if the hit is 100% identical to the query
        is_identical_seq = True 
        for i in range(dom_len):
            if pairwise_aln[i] != '.' and pairwise_aln[i] != dom_seq[i]:
                is_identical_seq = False
                break
        
        if is_identical_seq:
            continue
        
        if pairwise_aln not in out_alns:
            out_alns.append(pairwise_aln)
        
    ofile = open('bl_out2/'+dom_id,'w')
    ofile.write(dom_seq+'\n')
    for i in range(len(out_alns)):
        if i > 5000:
            break;
        ofile.write(out_alns[i]+'\n')
    ofile.close()

ifile.close()

    
    