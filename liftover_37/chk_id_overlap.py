#!/usr/bin/env python 
import sys

bim_file_name    = sys.argv[1]
strand_file_name = sys.argv[2]

# read SNP ids in the .bim file 
bfile=open(bim_file_name)
bim_snp_ids=[]
for line in bfile:
    bim_snp_ids.append(line.split()[1])
bfile.close()

# read SNP ids in the .strand file 
sfile=open(strand_file_name)
strand_snp_ids=[]
for line in sfile:
    strand_snp_ids.append(line.split()[0])
sfile.close()

# counting overlapping 
num_ids_overlapping    = 0
num_ids_only_in_bim    = 0
num_ids_only_in_strand = 0 

s_id_set=set(strand_snp_ids)
for snp_id in bim_snp_ids:
    if snp_id in s_id_set:
        num_ids_overlapping +=1
    if snp_id not in s_id_set:
        num_ids_only_in_bim +=1

b_id_set=set(bim_snp_ids)
for snp_id in strand_snp_ids:
    if snp_id not in b_id_set:
        num_ids_only_in_strand +=1

print len(bim_snp_ids), len(strand_snp_ids),\
      num_ids_only_in_bim, num_ids_overlapping, num_ids_only_in_strand
