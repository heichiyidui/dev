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

# get overlapping 
overlap_ids=[]
s_id_set=set(strand_snp_ids)
for snp_id in bim_snp_ids:
    if snp_id in s_id_set:
        overlap_ids.append(snp_id)

o_id_set=set(overlap_ids)

# read SNP allele in the strand file 
sfile=open(strand_file_name)
strand_allele={}
for line in sfile:
    snp_id=line.split()[0]
    if snp_id not in o_id_set:
        continue
    snp_allele=line.split()[5]
    strand_allele[snp_id]=snp_allele
sfile.close()

# read and check SNP allele in the bim file 
bfile=open(bim_file_name)
num_same_allele=0

for line in bfile:
    snp_id = line.split()[1]
    if snp_id not in o_id_set:
        continue
    a1=line.split()[4]+line.split()[5]
    a2=line.split()[5]+line.split()[4]
    
    if a1 == strand_allele[snp_id]:
        num_same_allele +=1
    if a2 == strand_allele[snp_id]:
        num_same_allele +=1
bfile.close()

print len(overlap_ids), num_same_allele

