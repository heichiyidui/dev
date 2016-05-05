################################################################################
#                Finding variants in PCSK9 affecting LDL-C                     #
################################################################################

################################################################################
# The PCSK9 gene blocks LDL (Low-density lipoprotein) receptors. Less LDL      #
# receptor on the surface of liver to remove LDL from bloodstream, leads to    #
# higher LDL cholesterol (LDL-C) concentrations.                               #
# PCSK9 is a fat-control drug target. It is well studied with caucasian        #
# populations.                                                                 #
# Now verify our CKB dataset on the detection of SNPs associating with LDL-C   #
# level in this gene.                                                          #
################################################################################


################################################################################
# 1. genotype data PCA

#######################################
# 1.1 the original set

# 32435 subjects in GWAS_sample_ascertainment.txt
# let's start from the stage3 set at
/kuser/shared/data/GWASphase12

tail -n +2 /kuser/shared/data/GWASphase12/stage3_mandatory_exclusions.txt | \
    awk '{print $1,$2,0,0,0,0}' > t.fam

plink --bfile /kuser/shared/data/GWASphase12/stage3 \
      --remove t.fam --make-bed --out ckb_ph12_s3
# 659231 variants and 32205 people
# 15599 '1' (male) and 16606 '2' (female)
# all with unknown father, mother and status
# no repeat in individual ids

# change the last column from '-9' to '0'

GWAS_sample_ascertainment.txt
# The list of sample ascertainments is from
# K:\kadoorie\Groups\Genetics\Data Archive\Project Sample Lists\Lists
# GWAS_SNPdata_samples.xlsx
# I replaced ' ' with '_' in the 'ascert.' and 'notes' columns.

# 10 RCs
# 1388 12 Qingdao
# 3310 16 Harbin
# 1173 26 Haikou
# 1705 36 Suzhou
# 2445 46 Liuzhou
# 4106 52 Sichuan
# 4703 58 Gansu
# 4341 68 Henan
# 3360 78 Zhejiang
# 5881 88 Hunan

#######################################
# 1.2 basic QC

tail -n +2 GWAS_sample_ascertainment.txt | awk '{print $1}' > t.ls
grab -f t.ls -c 2 ckb_ph12_s3.fam > t.fam

plink --bfile ckb_ph12_s3 --keep t.fam --make-bed --out ckb_ph12_s3_qc01
# 32109 individuals

# remove the bad SNPs found in manual clustering check
tail -n +2 /kuser/kuangl/dev/ckb_batch_check/manual_chk_res.table | \
    awk '{if ($2==0) print $1}' > to_rm_snp.ls
tail -n +2 /kuser/kuangl/dev/ckb_batch_check/plate_man_qc.table | \
    awk '{if ($2==0) print $1}' | awk -F"_" '{print $1}' >>  to_rm_snp.ls
# 14162 SNPs to be removed

plink --bfile ckb_ph12_s3_qc01 \
      --exclude to_rm_snp.ls   \
      --geno 0.05 \
      --maf  0.001 \
      --hwe 1e-8 midp \
      --make-bed --out ckb_ph12_s3_qc02

# 659231 SNPs, 32109 subjects in
#  13220 SNPs removed using the list 'to_rm_snp.ls'
#  49611 SNPs removed due to '--geno 0.05'
#  13099 SNPs removed due to '--hwe 1e-8 midp'
#  43128 SNPs removed due to '--maf  0.001'

# 540173 SNPs x 32109 people left

#######################################
# 1.3 LD-based pruning

plink  --bfile ckb_ph12_s3_qc02 \
    --autosome \
    --geno 0.01 \
    --hwe 1e-4 midp \
    --maf 0.05 \
    --indep-pairwise 1500 150 0.2

# 522460 SNPs loaded(autosome), 32109 subjects
#  15026 SNPs removed due to missing genotype data (--geno).
#   7298 SNPs removed due to Hardy-Weinberg exact test (--hwe).
# 164073 SNPs removed due to minor allele threshold(s)

# 336063 SNPs in, 216173 pruned, 119890 left

plink --bfile ckb_ph12_s3_qc02 \
    --extract plink.prune.in \
    --make-bed --out  ckb_ph12_pca
# 119890 SNPs, 32435 people

#######################################
# 1.4 split the plink files into 10 RC sets

for rc in 12 16 26 36 46 52 58 68 78 88 ; do
    tail -n +2 GWAS_sample_ascertainment.txt | \
        awk -v rc="$rc" '{if (substr($2,1,2)==rc) print $1}' > rc_id.ls
    grab -f rc_id.ls -c 2 ckb_ph12_pca.fam > t.fam
    plink --bfile ckb_ph12_pca \
          --keep t.fam \
          --make-bed --out ckb_ph12_pca_rc$rc
done

# Pretend all are control for the moment
for rc in 12 16 26 36 46 52 58 68 78 88 ; do
    awk '{$6=1; print $0}' ckb_ph12_pca_rc$rc.fam > t.fam
    mv t.fam  ckb_ph12_pca_rc$rc.fam
done

#######################################
# 1.5 families clusters


for rc in 12 16 26 36 46 52 58 68 78 88 ; do
    nohup plink --bfile ckb_ph12_pca_rc$rc \
          --genome --out rc$rc  &
done

# according to genome ibs distances,
# there are relatives in the cohorts, but no duplications now.

for rc in 12 16 26 36 46 52 58 68 78 88 ; do
    echo $rc
    tail -n +2 rc$rc.genome | awk '{if ($10>0.40) print $2 "\n" $4}' | \
        sort | uniq | wc
done

# some families to be removed during PCA
# using the threshold phat >= 0.125 for adding relatives

for ifile in *.genome ; do
    echo $ifile
    count_fam.py $ifile | sort -g | tail -n 3
done

# 12 Qingdao  1388      4 in the lastest family
# 16 Harbin   3310      3 in the lastest family
# 26 Haikou   1173     10 in the lastest family
# 36 Suzhou   1705    129 in the lastest family
# 46 Liuzhou  2445      4 in the lastest family
# 52 Sichuan  4106     12 in the lastest family
# 58 Gansu    4703     11 in the lastest family
# 68 Henan    4341      9 in the lastest family
# 78 Zhejiang 3360     18 in the lastest family
# 88 Hunan    5881     12 in the lastest family

# put the 129 ids of the Suzhou family into 'suzhou_fam.ls'



# for rc in 12 16 26 36 46 52 58 68 78 88 ; do
#     plink --bfile ckb_ph12_pca_rc$rc \
#           --read-genome rc$rc.genome \
#           --cluster \
#           --out rc$rc
# done
# all single clusters, mostly useless
# Should be some way to tweak it.

#######################################
# 1.6 First PCA

rc=36
~/bin/EIG/bin/smartpca.perl \
    -i ckb_ph12_pca_rc$rc.bed \
    -a ckb_ph12_pca_rc$rc.bim \
    -b ckb_ph12_pca_rc$rc.fam \
    -o t100_o.pca  \
    -p t100_o.plot \
    -e t100_o.eval \
    -l t100_o.log



grab -v -f suzhou_fam.ls -c 2  ckb_ph12_pca_rc36.fam | \
    awk '{print  $2, $5}' | \
    awk '{if ($2=="1") $2="M"; if ($2=="2") $2="F"; print $1,$2,"Control"}' | \
    sed 's/\ /\t/g'    > suzhou_no_fam.indiv

echo -e "suzhou_no_fam.indiv" > suzhou.pop

rc=36
~/bin/EIG/bin/smartpca.perl \
    -i ckb_ph12_pca_rc$rc.bed \
    -a ckb_ph12_pca_rc$rc.bim \
    -b ckb_ph12_pca_rc$rc.fam \
    -w suzhou_no_fam.indiv \
    -o t101_o.pca  \
    -p t101_o.plot \
    -e t101_o.eval \
    -l t101_o.log


# 601.57s user 65.66s system 180% cpu 6:10.50 total
# 6 min for 2000 individuals
# about 5 hours on 32435 individuals
# But out of memory on a PC with 12G memory.

# The GSL library is missing on the NC2! Emailed the administrator.
