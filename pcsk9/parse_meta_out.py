#!/usr/bin/env python3

# read SNP chr bp
snp_chr = {}
snp_bp  = {}
ifile=open('geno.bim')
for line in ifile:
    cols = line[:-1].split()
    SNP = cols[1]
    snp_chr[SNP] = cols[0]
    snp_bp [SNP] = cols[3]

ifile.close()

# read beta, se and p from linear results
snp_linear_beta = []
snp_linear_se   = []
snp_linear_p    = []

for i in range(5):
    snp_linear_p.append({})
    snp_linear_se.append({})
    snp_linear_beta.append({})

for i in range(5):
    ifile=open('st'+str(i+1)+'.raw.assoc.linear')
    ifile.readline()
    for line in ifile:
        cols = line[:-1].split()
        SNP = cols[1]
        snp_linear_beta[i][SNP]=cols[7]
        snp_linear_se  [i][SNP]=cols[8]
        snp_linear_p   [i][SNP]=cols[12]
    ifile.close()

# read meta beta, se, p and dir and Q
snp_meta_beta={}
snp_meta_se  ={}
snp_meta_p   ={}
snp_meta_q   ={}
snp_meta_dir ={}

ifile=open('pcsk9_direct1.tbl')
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    SNP=cols[0]
    snp_meta_beta[SNP] = cols[3]
    snp_meta_se  [SNP] = cols[4]
    snp_meta_p   [SNP] = cols[5]
    snp_meta_dir [SNP] = cols[6]
    snp_meta_q   [SNP] = cols[10]

ifile.close()

# print output
# title
print('SNP CHR BP',end = ' ')

print('str_1_p str_1_beta str_1_se',end=' ')
print('str_2_p str_2_beta str_2_se',end=' ')
print('str_3_p str_3_beta str_3_se',end=' ')
print('str_4_p str_4_beta str_4_se',end=' ')
print('str_5_p str_5_beta str_5_se',end=' ')

print('meta_beta meta_se meta_q meta_dir meta_p')

for SNP in snp_meta_p.keys():
    print(SNP,snp_chr[SNP],snp_bp[SNP],end=' ')
    for i in range(5):
        print(snp_linear_p   [i][SNP],  \
              snp_linear_beta[i][SNP],
              snp_linear_se  [i][SNP],end=' ')
    print(snp_meta_beta[SNP], \
          snp_meta_se  [SNP], \
          snp_meta_q   [SNP], \
          snp_meta_dir [SNP], \
          snp_meta_p   [SNP] )


