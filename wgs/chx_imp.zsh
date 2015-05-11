################################################################################
# notes and tutorials 
# http://www.nature.com/nrg/journal/v11/n7/box/nrg2796_BX3.html
# http://genome.sph.umich.edu/wiki/IMPUTE2:_1000_Genomes_Imputation_Cookbook
#
# http://mathgen.stats.ox.ac.uk/impute/example_one_phased_panel_chrX.html
# http://www.shapeit.fr/pages/m03_phasing/imputation.html

################################################################################
# 1. the original plink files 

# the merged UK set 

# in /scratch/home/ifoghbrc/INDEP_INTERNATIONAL_QC/UK_MERGE_QC 
# MERGE_UK_98_S.bim
# MERGE_UK_98_S.bed
# MERGE_UK_98_S.fam 

# need to apply QC for individuals
# UK-Cov.sample
# in /scratch/home/ifoghbrc/INTERNAT_SNPTEST/FAM_FILES

# need to apply QC for SNPs 
# --geno 0.01
# --hwe 0.000001

# select X, Y etc SNPS
awk '{if ($1>22) print $2}' \
 /scratch/home/ifoghbrc/INDEP_INTERNATIONAL_QC/UK_MERGE_QC/MERGE_UK_98_S.bim \
 > t_snp.ls 

# select individuals 
grep -v ID_1 \
  /scratch/home/ifoghbrc/INTERNAT_SNPTEST/FAM_FILES/MERGE_UK_AAO.sample \
  | grep -v "0 0 0 D" \
  | awk '{if (print $1,$2}' > t_id.ls 
  
# in the plink fam file, the cases all have the same individual id "1", 
# replace the individual ids with the different family ids 
awk '{$2=$1; print $0}' \
  /scratch/home/ifoghbrc/INDEP_INTERNATIONAL_QC/UK_MERGE_QC/MERGE_UK_98_S.fam \
  > t.fam

cp /scratch/home/ifoghbrc/INDEP_INTERNATIONAL_QC/UK_MERGE_QC/MERGE_UK_98_S.bim \
  t.bim 
cp /scratch/home/ifoghbrc/INDEP_INTERNATIONAL_QC/UK_MERGE_QC/MERGE_UK_98_S.bed \
  t.bed

plink --noweb --bfile t \
  --geno 0.01 \
  --hwe 0.000001 \
  --keep t_id.ls \
  --extract t_snp.ls \
  --make-bed --out uk

rm t.fam t.bed t.bim 

# other cohorts from 
# /home/ifoghbrc/Chr_X_QC

cp /home/ifoghbrc/Chr_X_QC/COR_chr.X_SAMPfltr_G_99.bim  co.bim 
cp /home/ifoghbrc/Chr_X_QC/D1_chr.X_SAMPfltr_G_99.bim   d1.bim 
cp /home/ifoghbrc/Chr_X_QC/D2_chr.X_SAMPfltr_G_99.bim   d2.bim 
cp /home/ifoghbrc/Chr_X_QC/IR_chr.X_SAMPfltr_G_99.bim   ir.bim 
cp /home/ifoghbrc/Chr_X_QC/M3_chr.X_SAMPfltr_G_99.bim   m3.bim 
cp /home/ifoghbrc/Chr_X_QC/SLA_chr.X_SAMPfltr_G_99.bim  sl.bim 

cp /home/ifoghbrc/Chr_X_QC/COR_chr.X_SAMPfltr_G_99.bed  co.bed 
cp /home/ifoghbrc/Chr_X_QC/D1_chr.X_SAMPfltr_G_99.bed   d1.bed 
cp /home/ifoghbrc/Chr_X_QC/D2_chr.X_SAMPfltr_G_99.bed   d2.bed 
cp /home/ifoghbrc/Chr_X_QC/IR_chr.X_SAMPfltr_G_99.bed   ir.bed 
cp /home/ifoghbrc/Chr_X_QC/M3_chr.X_SAMPfltr_G_99.bed   m3.bed 
cp /home/ifoghbrc/Chr_X_QC/SLA_chr.X_SAMPfltr_G_99.bed  sl.bed 

cp /home/ifoghbrc/Chr_X_QC/COR_chr.X_SAMPfltr_G_99.fam  co.fam 
cp /home/ifoghbrc/Chr_X_QC/D1_chr.X_SAMPfltr_G_99.fam   d1.fam 
cp /home/ifoghbrc/Chr_X_QC/D2_chr.X_SAMPfltr_G_99.fam   d2.fam 
cp /home/ifoghbrc/Chr_X_QC/IR_chr.X_SAMPfltr_G_99.fam   ir.fam 
cp /home/ifoghbrc/Chr_X_QC/M3_chr.X_SAMPfltr_G_99.fam   m3.fam 
cp /home/ifoghbrc/Chr_X_QC/SLA_chr.X_SAMPfltr_G_99.fam  sl.fam 

# 7 cohorts, co d1 d2 ir m3 sl uk
# 500 ~ 3500 individuals, 6000 ~ 13000 SNPs 
# coriell has no sex, use 
# /home/ifoghbrc/scratch/INDEP_INTERNATIONAL_QC/NIH_CORIELL/NIH_Corriel_TO_IMPUTE.fam
cp ~ifoghbrc/scratch/INDEP_INTERNATIONAL_QC/NIH_CORIELL/NIH_Corriel_TO_IMPUTE.fam \
  co.fam 

################################################################################
# 2. SNPs in the reference 
refdir=/scratch/data/1000genomes/ALL_1000G_phase1integrated_v3/\
ALL_1000G_phase1integrated_v3_impute

grep -v position $refdir/ALL_1000G_phase1integrated_v3_chrX_PAR1_impute.legend \
  | awk '{print $1}' > par1.ls 
grep -v position $refdir/ALL_1000G_phase1integrated_v3_chrX_PAR2_impute.legend \
  | awk '{print $1}' > par2.ls 
grep -v position $refdir/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend \
  | awk '{print $1}' > nonpar.ls 

#######################################
# grab is a simple python3 script in ~kuang/bin 
# 
# grab -f id.ls -c 2 in.table 
# means 
# first get the ids from the file "id.ls", 
# check the input table file "in.table",
# if for a line in "in.table", an id from "id.ls" is presented in the column 2,
# print the line to the standard ouput.

for cohort in co d1 d2 ir m3 sl uk; do 
    ~/bin/grab -f par1.ls -c 2 ${cohort}.bim | wc 
done 

for cohort in co d1 d2 ir m3 sl uk; do 
    ~/bin/grab -f par2.ls -c 2 ${cohort}.bim | wc 
done 

# We don't have any par1 or par2 SNPs.

# get the lists of nonpar chrx SNPs 
for cohort in co d1 d2 ir m3 sl uk; do 
    ~/bin/grab -f nonpar.ls -c 2 ${cohort}.bim |awk '{print $2}' >${cohort}_x.ls
done 
# In the 70432 SNPs in 7 cohorts, 69628 are nonpar chrx SNPs in the reference.

# extract the SNPs that also present in the reference
for cohort in co d1 d2 ir m3 sl uk; do 
    plink --noweb --bfile $cohort --extract ${cohort}_x.ls --make-bed \
    --out ${cohort}_01 
done 
# lost 804 SNPs, about 1%

# From chrx 2693518 to 154582606, about 152,000,000 bp. 
# We have around 10,000 SNPs on chr X, that's 1,000 per 15,000,000 bp 
# Less than what we need for imputation? 
# The reference has ~1,250,000 SNPs on chr X nonpar, around 100 times more.

# compared to chr 22 14438758 to 49565872, about 35,000,000 bp 
# about 9,000 SNPs on chr 22, about 1000 SNPs per 3,800,000 bp 

# will go on with imputation, expecting warnings of too few input SNPs from 
# the IMPUTE2 program. 

################################################################################
# 3. flip the strands, update to build37 positions 

cat ??_01.bim | awk '{print $5 "" $6}' |  sort | uniq -c

#    249 0A
#    274 0C
#    269 0G
#    248 0T
#   4269 AC
#  18962 AG
#      9 AT
#   3837 CA
#     16 CG
#   9594 CT
#  17760 GA
#     14 GC
#   2089 GT
#      6 TA
#   9692 TC
#   2340 TG

# 0A etc SNPs to be removed. 
# only a few AT TA CG and GC SNPs. remove them. 
for cohort in  co d1 d2 ir m3 sl uk; do 
    awk '{al = $5 "" $6; if ($5 == 0 || al == "AT" || al == "TA" || \
          al == "CG" || al == "GC") print $2}' ${cohort}_01.bim > t.ls 
    plink --noweb --bfile ${cohort}_01 \
        --exclude t.ls --make-bed --out ${cohort}_02
done 

cat ??_02.bim | awk '{print $2}' | sort | uniq > t.ls
wc t.ls 
# 14357 SNPs left from 14871 in the original bim files 

# download strand files from http://www.well.ox.ac.uk/~wrayner/strand/
# then 
cat BDCHP-1x10-HUMANHAP300v1-1_11219278_C-b37.strand    \
BDCHP-1X10-HUMANHAP550_11218540_C-b37.strand            \
Human610-Quadv1_B-b37-v2.strand                         \
Human660W-Quad_v1_A-b37.strand                          \
Human660W-Quad_v1_C-b37.strand                          \
HumanCNV370v1_C-b37-SourceStrand.strand                 \
HumanCNV370v1_C.strand                                  \
HumanHap300v2_A-b37.strand                              \
HumanHap550-2v3_B-b37.Illmn.strand                      \
HumanHap550-2v3_B-b37-SourceStrand.strand               \
HumanHap550-2v3_B-b37.strand                            > t.flip.strand

grab -f t.ls t.flip.strand > t_x.flip.strand
awk '{print $1}' t_x.flip.strand| sort | uniq | wc 
# for the 14357 SNPs, we have strand info for 14353

# use a simple python script to remove the redundency in the strand file 
# each SNP only presents once. 
t1.py 
#--------------------------------------
# #!/usr/bin/env python3 
# 
# ifile=open('t_x.flip.strand')
# snp_ids=set([])
# for line in ifile:
#     id = line.split()[0]
#     if id in snp_ids:
#         continue
#     snp_ids.add(id)
#     print(line[:-1])
# ifile.close()
#--------------------------------------

t1.py > x.strand 

wget http://www.well.ox.ac.uk/~wrayner/strand/update_build.sh
chmod +x update_build.sh

for cohort in  co d1 d2 ir m3 sl uk; do 
    update_build.sh ${cohort}_02 x.strand ${cohort}_03
done 

################################################################################
# 4. change to gen file format 

# 4.1 from binary plink to ped and map files
for cohort in  co d1 d2 ir m3 sl uk; do 
    plink --noweb --bfile ${cohort}_03 --recode --out ${cohort}_03
done

# 4.2 from plink ped and map to gen files 
for cohort in  co d1 d2 ir m3 sl uk; do 
    gtool -P --ped ${cohort}_03.ped --map ${cohort}_03.map \
              --og ${cohort}_03.gen --os ${cohort}_03.sample 
done 

# According to mathgen.stats.ox.ac.uk/impute/chromosome_X_file_formats.html
# Males have only two possible genotypes on chromosome X. The corresponding 
# genotypes in the gen files should be "1 0 0 " and "0 0 1 ".
# The heterozygous state should never be there. 
# Using a simple python script to check that.
t1.py 
#---------------------------------------
##!/usr/bin/env python3 
#import sys
#cohort=sys.argv[1]
#
#ifile=open(cohort+'_03.sample')
#gender=[]
#ifile.readline()
#ifile.readline()
#for line in ifile:
#    gender.append(line.split()[-2])
#ifile.close()
#
#ifile=open(cohort+'_03.gen')
#for line in ifile:
#    cols=line.split()
#    for i in range(len(gender)):
#        if gender[i]=='1':
#            print(cols[5+i*3 +1])
#ifile.close()
#
#---------------------------------------

for cohort in co d1 d2 ir m3 sl uk; do 
    t1.py $cohort | grep 1 | wc 
done 

#    289     289     578
#   6183    6183   12366
#    849     849    1698
#    232     232     464
#    380     380     760
#    943     943    1886
#    321     321     642

# So, most cohorts are fine. D1 is trouble here. 
# Check D1 males
# 2081870 "1 0 0 "
#    6183 "0 1 0 "
# 2148494 "0 0 1 "
#       1 "0 0 0 "
# Not too bad. About 0.15%. Leave it like that. 

################################################################################
# 5. imputation 
refdir=/scratch/data/1000genomes/ALL_1000G_phase1integrated_v3/\
ALL_1000G_phase1integrated_v3_impute

impute2 \
 -chrX \
 -m $refdir/genetic_map_chrX_nonPAR_combined_b37.txt \
 -h $refdir/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap \
 -l $refdir/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend \
 -g co_03.gen \
 -sample_g co_03.sample \
 -int 10.3e6 10.7e6 \
 -Ne 20000 \
 -o co_001.impute2

# from 2700157 to 148890965, 30 5x10^6 chunks

#--------------------------------------

qsub impu.sh co   1e6   6e6 001
qsub impu.sh co   6e6  11e6 002
qsub impu.sh co  11e6  16e6 003
qsub impu.sh co  16e6  21e6 004
qsub impu.sh co  21e6  26e6 005
qsub impu.sh co  26e6  31e6 006
qsub impu.sh co  31e6  36e6 007
qsub impu.sh co  36e6  41e6 008
qsub impu.sh co  41e6  46e6 009
qsub impu.sh co  46e6  51e6 010
qsub impu.sh co  51e6  56e6 011
qsub impu.sh co  56e6  61e6 012
qsub impu.sh co  61e6  66e6 013
qsub impu.sh co  66e6  71e6 014
qsub impu.sh co  71e6  76e6 015
qsub impu.sh co  76e6  81e6 016
qsub impu.sh co  81e6  86e6 017
qsub impu.sh co  86e6  91e6 018
qsub impu.sh co  91e6  96e6 019
qsub impu.sh co  96e6 101e6 020
qsub impu.sh co 101e6 106e6 021
qsub impu.sh co 106e6 111e6 022
qsub impu.sh co 111e6 116e6 023
qsub impu.sh co 116e6 121e6 024
qsub impu.sh co 121e6 126e6 025
qsub impu.sh co 126e6 131e6 026
qsub impu.sh co 131e6 136e6 027
qsub impu.sh co 136e6 141e6 028
qsub impu.sh co 141e6 146e6 029
qsub impu.sh co 146e6 151e6 030

# do the same for all 7 cohorts

# Jobs took around 1~3G memory. Asking for 5G should be safe. 
# Jobs run for around 2~3 hours. However, a few d2 jobs took about 24 hours. 

l -S impute1.e* | head 
l -S impute1.e* | tail
l -S *warnings  | head
l -S *warnings  | tail

# all empty files, no problem found 

for cohort in  co d1 d2 ir m3 sl uk; do 
    cat ${cohort}_???.impute2 > ${cohort}_04.gen 
    cp ${cohort}_03.sample ${cohort}_04.sample 
done 


awk '{print $2}' co_04.gen | sort | uniq -d 

# However, some positions of SNPs are wrong after updating to build37 
# e.g 
# rs1003590 48818427 vs 48818436
# rs1010973 101202727 vs 101202729

# end up with two entries for them in the impute2 files 
# --- rs1010973 101202727 A G 0 0 1 etc # imputed 
# 23 rs1010973 101202729 A G 0 0 1 etc # from original gen file with wrong pos 

# note the imputed SNP entries start with '---'
# keep the original SNP, remove the imputed one, 
# and update the position of the original one. 

# then change all '---' to '23'
t1.py 
#--------------------------------------
##!/usr/bin/env python3 
#import sys
#
#cohort=sys.argv[1]
#
## read build37 SNP positions
#ifile=open('/scratch/data/1000genomes/ALL_1000G_phase1integrated_v3/\
#ALL_1000G_phase1integrated_v3_impute/\
#ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend')
#
#ifile.readline()
#snp_pos={}
#for line in ifile:
#    cols=line.split()
#    snp_pos[cols[0]]=cols[1]
#ifile.close()
#
## get the list of original SNPs 
#org_snp=set([])
#ifile=open(cohort+'_03.gen')
#for line in ifile:
#    cols=line.split()
#    org_snp.add(cols[1])
#ifile.close()
#
## read and write gen files
#ifile=open(cohort+'_04.gen')
#for line in ifile:
#    cols=line.split()
#    if cols[0] == '---':
#        if cols[1] in org_snp:
#            continue; # ignore the imputed genotypes when the original is there 
#        cols[0] = '23'
#        print(' '.join(cols))
#    else: 
#        cols[2] = snp_pos[cols[1]]
#        print(' '.join(cols))
#ifile.close()
#--------------------------------------

# rs145285390 is alwasys the trouble. 
# it is located at 116000000. There are two impute jobs covering it. 
# qsub impu.sh co 111e6 116e6 023
# qsub impu.sh co 116e6 121e6 024
# It will be written twice, to two .impute2 files 
# Need to remove one entry of it from the gen files. 
# The output gen files are always of 1223763 lines. 
# The SNP rs145285390 is always located at line 927811 and 927812

for cohort in  co d1 d2 ir m3 sl uk; do 
    t1.py $cohort > t.out
    head -n 927811 t.out > $cohort\_05.gen 
    tail -n 295951 t.out >> $cohort\_05.gen 
    cp ${cohort}_04.sample ${cohort}_05.sample 
done 





