################################################################################
#                Finding variants in PCSK9 affecting LDL-C                     #
################################################################################

################################################################################
#                                                                              #
# The PCSK9 gene blocks LDL (Low-density lipoprotein) receptors. Less LDL      #
# receptor on the surface of liver to remove LDL from bloodstream, leads to    #
# higher LDL cholesterol (LDL-C) concentrations.                               #
# PCSK9 is a fat-control drug target. It is well studied with caucasian        #
# populations.                                                                 #
#                                                                              #
# Now verify our CKB dataset on the detection of SNPs associating with LDL-C   #
# level in this gene.                                                          #
#                                                                              #
################################################################################

# K:\kadoorie\Groups\Genetics\PROJECTS\PCSK9

################################################################################
# 1. genotype data PCA

#######################################
# 1.1 the original set

# 32435 subjects in GWAS_sample_ascertainment.txt from
# K:\kadoorie\Groups\Genetics\Data Archive\Project Sample Lists\Lists
# I replaced ' ' with '_' in the 'ascert.' and 'notes' columns.

# let's start from the stage3 set at
/kuser/shared/data/GWASphase12

plink --bfile /kuser/shared/data/GWASphase12/stage3 \
      --remove /kuser/shared/data/GWASphase12/stage3_mandatory_exclusions.txt \
      --make-bed --out ckb_ph12_s3

# 659231 variants and 32205 people
# all with unknown father, mother and status
# no repeat in individual ids

# make-sure all subjects have ascertainments
tail -n +2 GWAS_sample_ascertainment.txt | awk '{print $1}' > t.ls
grab -v -f t.ls -c 2 ckb_ph12_s3.fam > t.fam
# 96 subjects to be removed

plink --bfile ckb_ph12_s3 \
      --remove t.fam \
      --autosome \
      --geno 0.05 \
      --maf  0.001 \
      --make-bed --out ckb_ph12_s3_qc01
# 542534 variants and 32109 people
# 51983 variants removed due to missing genotype data (--geno).
# 42153 variants removed due to minor allele threshold(s)

awk '{print $5}' ckb_ph12_s3_qc01.fam | sort | uniq -c
# 15599 '1' (male) and 16606 '2' (female)


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



# skip the PCA and QC part below. Go section 3.




# remove the bad SNPs found in manual clustering check
tail -n +2 /kuser/kuangl/dev/ckb_batch_check/manual_chk_res.table | \
    awk '{if ($2==0) print $1}' > to_rm_snp.ls
tail -n +2 /kuser/kuangl/dev/ckb_batch_check/plate_man_qc.table | \
    awk '{if ($2==0) print $1}' | awk -F"_" '{print $1}' >>  to_rm_snp.ls
# 14162 SNPs to be removed

plink --bfile ckb_ph12_s3_qc01 \
      --exclude to_rm_snp.ls   \
      --autosome \
      --geno 0.05 \
      --maf  0.001 \
      --hwe 1e-8 midp \
      --make-bed --out ckb_ph12_s3_qc02

# 636670 SNPs loaded (autosome), 32109 subjects
#  12346 SNPs removed using the list 'to_rm_snp.ls'
#  47053 SNPs removed due to '--geno 0.05'
#  12947 SNPs removed due to '--hwe 1e-8 midp'
#  41478 SNPs removed due to '--maf  0.001'

# 522846 SNPs x 32109 people left

#######################################
# 1.3 LD-based pruning

plink  --bfile ckb_ph12_s3_qc02 \
    --geno 0.01 \
    --hwe 1e-4 midp \
    --maf 0.05 \
    --indep-pairwise 1500 150 0.2

# 522846 SNPs loaded, 32109 subjects
#  15412 SNPs removed due to missing genotype data (--geno).
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
# edit the last column of ckb_ph12_pca_rc36.fam
# so that we have 'fam_1' or 'no_rel' in the last column

echo "no_rel" > t.pop
# population list file for smartpca on rc36

# for rc in 12 16 26 36 46 52 58 68 78 88 ; do
#     plink --bfile ckb_ph12_pca_rc$rc \
#           --read-genome rc$rc.genome \
#           --cluster \
#           --out rc$rc
# done
# all single clusters, mostly useless
# Should be some way to tweak it.

#######################################
# 1.6 first round region PCA

# no rc36
for rc in 12 16 26 46 52 58 68 78 88 ; do
    ~/bin/EIG/bin/smartpca.perl \
        -i ckb_ph12_pca_rc$rc.bed \
        -a ckb_ph12_pca_rc$rc.bim \
        -b ckb_ph12_pca_rc$rc.fam \
        -o rc$rc.pca  \
        -p rc$rc.plot \
        -e rc$rc.eval \
        -l rc$rc.log
done

for rc in 36 ; do
    ~/bin/EIG/bin/smartpca.perl \
        -i ckb_ph12_pca_rc$rc.bed \
        -a ckb_ph12_pca_rc$rc.bim \
        -b ckb_ph12_pca_rc$rc.fam \
        -w t.pop \
        -o rc$rc.pca  \
        -p rc$rc.plot \
        -e rc$rc.eval \
        -l rc$rc.log
done

# 6 min for 2000 individuals
# about 5 hours on 32435 individuals
# But out of memory on a PC with 12G memory.

# The GSL library is missing on the NC2! Emailed the administrator.

# So far only rc 46 PCA shows some internal structure there.
# Use the first 10 global PCs, from stage4_pca.xlsx

ckb_global_pca.txt

################################################################################
# 2 LDL-C measurements

# from
# file://K:\kadoorie\Groups\CKB-Statistics\Biochemistry Data\20160419_final_data
# 20160426_101958_biochemistry_result_grid.xls

# about 20K subjects
# removed the 'diluted'
# removed the 'Result may be affected by lipemia' etc

pheno.ods
# we need their DOB.
# we need to add sttroke and statin use.

# checking different covariants
plot_pheno.R

# put the ldl phenotype into ldl.pheno

# rc cov as 9 binary factors
# rc_map = {'12':0, '16':1, '26':2, '36':3, '46':4,
#           '52':5, '58':6, '68':7, '78':8}

# ascertainment cov as 4 binary factors
# rc_map = {'control':0, 'ICH':1, 'IS':2, 'MI_IHD':3}

# put the covariants into ldl.cov

################################################################################
# 3 genotype file

# almost without any QC

plink --bfile ckb_ph12_s3_qc01 \
      --chr 1 --from-mb 55.44 --to-mb 56.01 \
      --keep ldl.pheno \
      --make-bed --out geno
# 159 SNPs, 17458 subjects


################################################################################
# 4. assoc

plink --bfile geno \
      --pheno ldl.pheno \
      --linear --ci 0.95

# 1   AX-83389438   55509585    T        ADD
# 17455   -0.09347 0.007536  -0.1082 -0.07869        -12.4    3.588e-35

# with all covariants
plink --bfile geno \
      --pheno ldl.pheno \
      --covar ldl.cov \
      --linear --ci 0.95

grep ADD plink.assoc.linear | sort -g -k 12 | head -n 22
head  -n 1 plink.assoc.linear > t.in
grep ADD plink.assoc.linear | sort -g -k 12 | grep -v NA  >> t.in
# then plot it

#   1   AX-83389438   55509585    T        ADD
# 17455      -0.38  0.03056  -0.4399  -0.3201       -12.43    2.419e-35
# MAF = 0.013


# check ld
plink --bfile geno \
       --r2 --ld-snp AX-39912161 \
       --ld-window-r2 0 \
       --ld-window 99999 \
       --ld-window-kb 77000
# nothing is in LD with AX-83389438

sort -g -k 7 plink.ld
#   1     55509585   AX-83389438      1     55601339   AX-11187195    0.0140177
# highest of all but itself is 0.0140177

plink --bfile geno \
       --r2 --ld-snp AX-11576926 \
       --ld-window-r2 0 \
       --ld-window 99999 \
       --ld-window-kb 77000
# not too bad

plink --bfile geno \
       --r2 --ld-snp AX-39912161 \
       --ld-window-r2 0 \
       --ld-window 99999 \
       --ld-window-kb 77000
#############
# chk the clustering
