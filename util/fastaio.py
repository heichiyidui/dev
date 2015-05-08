#!/usr/bin/env python3 

def read_fasta(ifile_name):
    ''' To read FASTA entries from an input file. 
    A list of sequences will be returned. 
    Each sequence should be two strings combined by a endl,
    the first being the title line, the second the amino acid sequence.'''
    
    ifile=open(ifile_name)
    lines=ifile.read().split('\n')
    ifile.close()
    
    seqs=[]
    seq_buffer=''
    
    for line in lines:
        if line.strip()=='':
            continue
        if line.startswith('>'): # found another title line 
            if seq_buffer != '' :
                seqs.append(seq_buffer)
            seq_buffer=line+'\n'
            continue
        seq_buffer+=''.join(line.split(' ')).upper()
        # to remove spaces in the seq lines
        # and use uppercase only
    
    if seq_buffer != '':
        seqs.append(seq_buffer)
    
    return seqs

def write_fasta(ofile_name,seqs):
    ''' To write FASTA entries to an output file.
    A list of sequences (for format see read_fasta())
    will be written to the output file with file name ofile_name'''
    
    ofile=open(ofile_name,'w')
    for seq in seqs:
        (seq_title,seq_seq)=seq.split('\n')
        ofile.write(seq_title+'\n')
        for i in range(0,len(seq_seq),80):
            ofile.write(seq_seq[i:i+80]+'\n')
    ofile.close()

if __name__ == '__main__':
    seqs=read_fasta('caln/102mA00')
    for seq in seqs:
        print(seq)
    write_fasta('t.out',seqs[1:])

