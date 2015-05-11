#!/usr/bin/env python3

import sys

PIECES_SNP_SIZE  = 1200     # roughly that many SNPs in any pieces
PIECES_LENGTH    = 4500000  # rough piece size in bp
BUFFER_SIZE      = 250000   # default buffer size
MIN_PIECE_SIZE   = 50       # in a chunk, should have at least 51 SNPs 
#######################################
# to read SNP positions and chr 

ifile=open(sys.argv[1]) # input bim or map file
#ifile=open('pe19_22.bim')
snp_pos=[]  # SNP bp position

line=ifile.readline()
cols=line.split()
chr='{:02d}'.format(int(cols[0])) # 01 02 ... 22

snp_pos.append(int(cols[3]))

for line in ifile:
    snp_pos.append(int(line.split()[3]))
ifile.close()

#######################################
# to cut chromsomes into chunks
# SNPs in a chunk can not be seperated by more than PIECES_LENGTH / 4

chunks=[]
chunk=[snp_pos[0]]
for snp_po in snp_pos[1:]:
    if snp_po - chunk[-1] > PIECES_LENGTH / 4:
        chunks.append(chunk)
        chunk=[snp_po]
    else:
        chunk.append(snp_po)
chunks.append(chunk)

# print('\nchr:',chr, '\tnum_of_chunks:',len(chunks))

#######################################
# to cut chunks into pieces

pieces=[]
for chunk in chunks:
    # print('\tchunk_begin:',chunk[0],\
    #       '\tchunk_end:',  chunk[-1],\
    #       '\tchunk_length:',chunk[-1]-chunk[0]+1,\
    #       '\tSNPs_in_chunk:',len(chunk))
    if len(chunk) < MIN_PIECE_SIZE:
        continue
    len_chunk=chunk[-1]-chunk[0]
    num_pieces=len_chunk//PIECES_LENGTH + 1
    # print('\t\tnum_of_pieces_in_chunk:',num_pieces,\
    #       '\taverage_SNPs_in_pieces:',len(chunk)//num_pieces)
    
    size_pieces=len(chunk)//num_pieces
    for i in range(0,len(chunk),size_pieces):
        pieces.append(chunk[i:i+size_pieces])

for i in range(7):
    new_pieces=[]
    for piece in pieces:
        if len(piece) < MIN_PIECE_SIZE:
            continue
        if piece[-1] - piece[0] < PIECES_LENGTH:
            new_pieces.append(piece)
            continue
        
        dist=[] # distances between snps
        for j in range(len(piece)-1):
            dist.append(piece[j+1]-piece[j])
        sep=dist.index(max(dist))# the seperation point
        
        if sep < len(piece)*0.1 or sep > len(piece)*0.9:
            sep=len(piece)//2
        
        new_pieces.append(piece[0:sep+1])
        new_pieces.append(piece[sep+2:])
    pieces=new_pieces

#######################################
# to get the begins and ends of pieces

begins=[]
ends=[]
begins.append(max(pieces[0][0]-BUFFER_SIZE,1))

for i in range(len(pieces)-1):
    gap_begin = pieces[i][-1]
    gap_end   = pieces[i+1][0]
    
    ends.append  (min( (gap_begin+gap_end)//2  ,gap_begin+BUFFER_SIZE))
    begins.append(max( (gap_begin+gap_end)//2+1,gap_end  -BUFFER_SIZE))
ends.append(pieces[-1][-1]+BUFFER_SIZE)

#######################################
# to output the command

for i in range(len(begins)):
    piece_code='{:03d}'.format(i)
    print('qsub impu.sh',chr,piece_code,begins[i],ends[i])
