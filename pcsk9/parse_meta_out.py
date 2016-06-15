#!/usr/bin/env python3

#######################################
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

#######################################
# read SNP MAF
snp_maf ={}
ifile=open('plink.frq')
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    snp_maf[cols[1]] = cols[4]
ifile.close()

#######################################
# read beta, se and p from linear results
snp_linear_beta = []
snp_linear_se   = []
snp_linear_p    = []

for i in range(6):
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

# stratum 6 is a rint set
i=5
ifile=open('st'+str(i+1)+'.rint.assoc.linear')
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    SNP = cols[1]
    snp_linear_beta[i][SNP]=cols[7]
    snp_linear_se  [i][SNP]=cols[8]
    snp_linear_p   [i][SNP]=cols[12]
ifile.close()

#######################################
# read meta a1, a2, beta, se, p and dir and Q
snp_meta_0_a1  ={}
snp_meta_0_a2  ={}
snp_meta_0_beta={}
snp_meta_0_se  ={}
snp_meta_0_p   ={}
snp_meta_0_q   ={}
snp_meta_0_dir ={}

ifile=open('pcsk9_direct1.tbl')
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    SNP=cols[0]
    snp_meta_0_a1  [SNP] = cols[1].upper()
    snp_meta_0_a2  [SNP] = cols[2].upper()
    snp_meta_0_beta[SNP] = cols[3]
    snp_meta_0_se  [SNP] = cols[4]
    snp_meta_0_p   [SNP] = cols[5]
    snp_meta_0_dir [SNP] = cols[6]
    snp_meta_0_q   [SNP] = cols[10]
ifile.close()

snp_meta_1_a1  ={}
snp_meta_1_a2  ={}
snp_meta_1_beta={}
snp_meta_1_se  ={}
snp_meta_1_p   ={}
snp_meta_1_q   ={}
snp_meta_1_dir ={}

ifile=open('pcsk9_all1.tbl')
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    SNP=cols[0]
    snp_meta_1_a1  [SNP] = cols[1].upper()
    snp_meta_1_a2  [SNP] = cols[2].upper()
    snp_meta_1_beta[SNP] = cols[3]
    snp_meta_1_se  [SNP] = cols[4]
    snp_meta_1_p   [SNP] = cols[5]
    snp_meta_1_dir [SNP] = cols[6]
    snp_meta_1_q   [SNP] = cols[10]
ifile.close()


#######################################
# print output
# title
print('SNP CHR BP MAF',end = ' ')

# print('str_1_p str_1_beta str_1_se',end=' ')
# print('str_2_p str_2_beta str_2_se',end=' ')
# print('str_3_p str_3_beta str_3_se',end=' ')
# print('str_4_p str_4_beta str_4_se',end=' ')
# print('str_5_p str_5_beta str_5_se',end=' ')

print('beta_0 se_0 q_0 dir_0 p_0', end = ' ')
print('beta_1 se_1 q_1 dir_1 p_1')

for SNP in snp_meta_0_p.keys():
    if SNP not in snp_meta_1_p.keys():
        continue
    print(SNP,snp_chr[SNP],snp_bp[SNP],snp_maf[SNP],end=' ')

    # for i in range(5):
    #     print(snp_linear_p   [i][SNP],  \
    #           snp_linear_beta[i][SNP],
    #           snp_linear_se  [i][SNP],end=' ')
    print(snp_meta_0_beta[SNP], \
          snp_meta_0_se  [SNP], \
          snp_meta_0_q   [SNP], \
          snp_meta_0_dir [SNP], \
          snp_meta_0_p   [SNP], end =' ' )
    print(snp_meta_1_beta[SNP], \
          snp_meta_1_se  [SNP], \
          snp_meta_1_q   [SNP], \
          snp_meta_1_dir [SNP], \
          snp_meta_1_p   [SNP] )


