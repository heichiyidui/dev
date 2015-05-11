################################################################################
# polygenic score using the GCTA package                                       #
################################################################################

################################################################################
# 1. get the plink binary files                                                #
################################################################################

cd /home/klinbrc/scratch/dev/gcta # working directory

peddir=~ifoghbrc/scratch/CGTA-POLY/PLINK-FORMAT 
# the directory of orginal ped and map files

#######################################
# to get the tped and tfam files
for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for chr in 01 02 03 04 05 06 07 08 09 10 \
               11 12 13 14 15 16 17 18 19 20 21 22; do 
        plink --noweb --file ${peddir}/${cohort}_IMP_chr_${chr} 
              --recode --transpose --out ${cohort}_chr_${chr}
    done
done
# That's 8x22 plink commands. They can and should be splitted into 8x22
# queue jobs

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    cat ${cohort}_chr_??.tped > ${cohort}.tped
    cp ${cohort}_chr_22.tfam ${cohort}.tfam
done
# 8 slow 'cat' commands, can be splited to 8 queue jobs

#######################################
# to remove duplications in SNP ids and bp positions from the tped files

# rmSNP_POS_dup.py:

#--------------------------------------
#!/usr/bin/env python2.7
import sys
 
ifile=open(sys.argv[1])
ids=set([])
bps=set([])
for line in ifile:
    cols=line.split()
    (id,bp)=(cols[1],cols[0]+'_'+cols[3]) # bp is actually chr_bp here
    if id in ids:
        continue;
    if bp in bps:
        continue 
    print(line[:-1])
    ids.add(id)
    bps.add(bp)
ifile.close()
#--------------------------------------

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    ./rmSNP_POS_dup.py ${cohort}.tped > temp_${cohort}.tped; 
    mv temp_${cohort}.tped ${cohort}.tped
done
# again, should be 8 queue jobs

#######################################
# to get the plink binary files

mkdir bed

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    plink --noweb --tfile ${cohort} --make-bed --out bed/${cohort}
done
# again, should be 8 queue jobs
# plink asks for a lot of memory for the big sets

# the cohorts have 420 (IR) to 3536 (SLA) individuals
# the cohorts have 3.4 (M3) to 4.8 (IR) million SNPs

#######################################
# clean up the tped, tfam and plink log files...
rm *.tped 
rm *.tfam
rm *.log

# wc bed/*.fam
#   507   3042  12617 bed/COR.fam
#   852   5112  25600 bed/D1.fam
#  2599  15594 100629 bed/D2.fam
#   420   2520  10287 bed/IR.fam
#  2449  14694  71178 bed/M3.fam
#   653   3918  11744 bed/NIT.fam
#  3536  21216  79242 bed/SLA.fam
#  2209  13254  40832 bed/UK.fam
# 13225  79350 352129 total

# wc bed/*.bim
#  4594687  27568122 129291890 bed/COR.bim
#  3968571  23811426 111637724 bed/D1.bim
#  4067660  24405960 114438778 bed/D2.bim
#  4824119  28944714 135760869 bed/IR.bim
#  3448792  20692752  96988661 bed/M3.bim
#  4361134  26166804 122667534 bed/NIT.bim
#  4549968  27299808 128008282 bed/SLA.bim
#  4657408  27944448 131082315 bed/UK.bim
# 34472339 206834034 969876053 total


################################################################################
# 2. get the grm files with GCTA                                               #
################################################################################


#######################################
# get grm files 

#  download the current version of gcta
wget http://www.complextraitgenomics.com/software/gcta/gcta_1.04.zip
unzip gcta_1.04.zip 
# install the executable
mv gcta64 ~/bin
# cleaning up 
rm gcta_1.04.zip gcta_mac COPYING.txt test.bed test.bim test.fam \
    README.txt test.phen gcta gcta.exe

mkdir grm

# t.sh: 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m e
#$ -N gcta
#$ -cwd
#$ -q long.q,bignode.q,short.q

~/bin/gcta64 --bfile bed/$1 \
    --autosome --make-grm  \
    --save-ram --grm-cutoff 0.025 \
    --out grm/$1
#--------------------------------------

# A problem of memory usage
# gcta uses a LOT of memory
# There are 3536 individuals, 4.5 million SNPs in the full SLA set. 
# gcta used 64G memory for the calculation of grm 

# Assuming the memory usage is linear to the product of individual and SNP 
# numbers, we don't need to worry about the fist four cohorts

for cohort in D1 COR IR NIT ; do 
    qsub -l h_vmem=20G t.sh ${cohort} 
done

for cohort in SLA D2 M3 UK ; do 
    qsub -l h_vmem=70G t.sh ${cohort} 
done

ll -S gcta.e* | head

################################################################################
# 3. Estimation of the variance explained by the SNPs                          #
################################################################################

#PCs

#SLAGEN  4  3536
#D1      0  852
#D2      4  2599
#CORRIEL 0  507
#IRISH   1  420
#MGH     2  2449
#NIT     1  653
#UK      1  2209

#######################################
# get phenotype files 
rm t.ls
for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    awk '{print $6}' bed/${cohort}.fam >> t.ls
done
sort t.ls | uniq 
# only 1 and -9 in the phenotype columns of the family files

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    awk '{if ($6==-9) $6=0; print $1,$2,$6}' bed/${cohort}.fam \
        > ${cohort}_cc.phen
done

#######################################
# get the ratio of genetic variance to phenotypic variance

mkdir hsq

# two cohorts with no pc correction

for cohort in D1 COR; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.0001 \
           --out hsq/${cohort}_${p_thres}
    done
done

# for cohorts with PC corrections
# Why ID_1 and ID_2 are in different orders in .fam and sample files?
mkdir sample
mv *.sample sample

# 4 PCs
for cohort in SLA D2; do 
    grep -v "ID_1 ID_2" sample/${cohort}*.sample | grep -v "0 0 0" | \
        awk '{print $2,$1,$5,$6,$7,$8}' > ${cohort}.pc
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.0001 \
           --out hsq/${cohort}_${p_thres}
    done
done

# 2 PCs
for cohort in M3 ; do
    grep -v "ID_1 ID_2" sample/${cohort}*.sample | grep -v "0 0 0" | \
        awk '{print $2,$1,$5,$6}' > ${cohort}.pc
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.0001 \
           --out hsq/${cohort}_${p_thres}
    done
done

# 1 PC
for cohort in IR NIT UK ; do
    grep -v "ID_1 ID_2" sample/${cohort}*.sample | grep -v "0 0 0" | \
        awk '{print $2,$1,$5}' > ${cohort}.pc
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.0001 \
           --out hsq/${cohort}_${p_thres}
    done
done

cd hsq 
grep Vp_L *.hsq | sed 's/.hsq:V(G)\/Vp_L//'  | awk '{print $1 "\t" $2}'
cd ..

# see vp_l.xls for the results


################################################################################
# 4. do it without chr 9 and try --prevalence 0.00001                          #
################################################################################

#######################################
# 4.1 to get the plink binary files 
mkdir bed_no9

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        awk '{if ($1==9) print $2}' bed/${cohort}_${p_thres}.bim \
            > t_9_${cohort}_${p_thres}.ls
        plink --noweb --bfile bed/${cohort}_${p_thres} \
            --exclude t_9_${cohort}_${p_thres}.ls \
            --make-bed --out bed_no9/${cohort}_${p_thres}
        rm t_9_${cohort}_${p_thres}.ls
    done
done

# OR to put the command lines into t.sh with cohort=$1 p_thres=$2
# and run it like 
for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        qsub -l h_vmem=5G t.sh ${cohort} ${p_thres}
    done
done

rm bed_no9/*.log bed_no9/*.nosex

mkdir grm_no9

# put the command line into t.sh 
~/bin/gcta64 --bfile bed_no9/$1 --autosome --save-ram --make-grm \
 --out grm_no9/$1

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for p_thres in 0.001 0.002 0.005 0.01 ; do
        qsub t.sh ${cohort}_${p_thres};
    done
done

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for p_thres in 0.05; do
        qsub -l h_vmem=5G t.sh ${cohort}_${p_thres};
    done
done

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for p_thres in 0.1; do
        qsub -l h_vmem=10G t.sh ${cohort}_${p_thres};
    done
done

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for p_thres in 0.2; do
        qsub -l h_vmem=20G t.sh ${cohort}_${p_thres};
    done
done

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for p_thres in 0.3; do
        qsub -l h_vmem=26G t.sh ${cohort}_${p_thres};
    done
done

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for p_thres in 0.4 0.5 ; do
        qsub -l h_vmem=50G t.sh ${cohort}_${p_thres};
    done
done

for cohort in D1 COR IR NIT ; do 
    for p_thres in 1 ; do
        qsub -l h_vmem=20G t.sh ${cohort}_${p_thres};
    done
done

for cohort in SLA D2 M3 UK ; do 
    for p_thres in 1 ; do
        qsub -l h_vmem=67G t.sh ${cohort}_${p_thres};
    done
done

mkdir hsq_no9_45

# two cohorts with no pc correction
for cohort in D1 COR; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm_no9/${cohort}_${p_thres} \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00005 \
           --out hsq_no9_45/${cohort}_${p_thres}
    done
done

# 4 PCs
for cohort in SLA D2; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm_no9/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00005 \
           --out hsq_no9_45/${cohort}_${p_thres}
    done
done

# 2 PCs
for cohort in M3 ; do
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm_no9/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00005 \
           --out hsq_no9_45/${cohort}_${p_thres}
    done
done

# 1 PC
for cohort in IR NIT UK ; do
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm_no9/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00005 \
           --out hsq_no9_45/${cohort}_${p_thres}
    done
done


#######################################
# 4.2 prevalence 0.00001

mkdir hsq_5

# two cohorts with no pc correction
for cohort in D1 COR; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00001 \
           --out hsq_5/${cohort}_${p_thres}
    done
done

# 4 PCs
for cohort in SLA D2; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00001 \
           --out hsq_5/${cohort}_${p_thres}
    done
done

# 2 PCs
for cohort in M3 ; do
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00001 \
           --out hsq_5/${cohort}_${p_thres}
    done
done

# 1 PC
for cohort in IR NIT UK ; do
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00001 \
           --out hsq_5/${cohort}_${p_thres}
    done
done

cd hsq_5 
grep Vp_L *.hsq | sed 's/.hsq:V(G)\/Vp_L//'  | awk '{print $1 "\t" $2}'
cd ..

#######################################
# 4.3 prevalence 0.00005

mkdir hsq_45

# two cohorts with no pc correction
for cohort in D1 COR; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00005 \
           --out hsq_45/${cohort}_${p_thres}
    done
done

# 4 PCs
for cohort in SLA D2; do 
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00005 \
           --out hsq_45/${cohort}_${p_thres}
    done
done

# 2 PCs
for cohort in M3 ; do
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00005 \
           --out hsq_45/${cohort}_${p_thres}
    done
done

# 1 PC
for cohort in IR NIT UK ; do
    for p_thres in 0.001 0.002 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 1 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence 0.00005 \
           --out hsq_45/${cohort}_${p_thres}
    done
done

cd hsq_45
grep Vp_L *.hsq | sed 's/.hsq:V(G)\/Vp_L//'  | awk '{print $1 "\t" $2}'
cd ..

###############################################################################
# 5. pi_hat clustering of individuals.                                        #
###############################################################################

grep "__" bed/*_0.001.fam
# "__" is not presented in the individual ids 
awk '{print $1 "__" $2}' bed/*_0.001.fam  | sort | uniq -d
# no duplications in the individual ids

grep -v eig merge_8_LD_PIHAT_Clean.evec \
 | awk '{print $1,$2,$3,$4}' | sed 's/:/__/' > t.evec

awk '{print $2}' \
     ~ifoghbrc/scratch/PLINK_QC_COMBINED_COHORT/merge_8_LD.bim > t.ls

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    echo ${cohort}
    awk '{print $1 "__" $2}' bed/${cohort}.fam  > ind.ls 
    grab -f ind.ls t.evec | wc
done
# lots of troubles

#######################################
# 5.1 to get the combined LD-free set 
awk '{print $2}' \
   ~ifoghbrc/scratch/PLINK_QC_COMBINED_COHORT/merge_8_LD.bim > t.ls
# 70647 LD-free SNPs

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    echo ${cohort}
    plink --noweb --bfile bed/${cohort} --extract t.ls \
        --make-bed --out bed/${cohort}_ld
done
for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    qsub -l h_vmem=6G t2.sh ${cohort}
done

awk '{print $2}' bed/*_ld.bim | sort | uniq -c \
  | awk '{if ($1==8) print $2}' > t2.ls
wc t2.ls
# 64316 LD-free SNPs presented in all 8 QCed cohorts

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    echo ${cohort}
    plink --noweb --bfile bed/${cohort}_ld --extract t2.ls \
        --make-bed --out bed/${cohort}_ld
done

for cohort in D1 D2 IR M3 NIT SLA UK; do 
    awk '{print $1,$2,$3,$4}' bed/COR_ld.bim > t1.map
    awk '{print $1,$2,$3,$4}' bed/${cohort}_ld.bim > t2.map
    diff t1.map t2.map | wc
done
# all 0, the LD-free sets now have identical maps

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    plink --noweb --bfile bed/${cohort}_ld --recode \
        --out ${cohort}_ld
done

cat COR_ld.ped D1_ld.ped D2_ld.ped IR_ld.ped M3_ld.ped \
    NIT_ld.ped SLA_ld.ped UK_ld.ped > m8.ped 
cp COR_ld.map m8.map

plink --noweb --file m8 --make-bed --out bed/m8
# 64316 LD-free SNPs, 13225 individuals

awk '{ if ($6==1) $6=2; if ($6==-9) $6=1; print $0}' bed/m8.fam > t.out
mv t.out bed/m8.fam 
# correct the phenotypes

#######################################
# 5.2 eigenstrat

mkdir pca
cd pca
mv ../bed/m8.??? .

R --vanilla --slave --args stem=m8 nsnpldregress=5 < EIGENSOFTplus_v12.txt


qsub t2.sh -l h_vmem=32G
# For whatever reason, it keeps crash, giving a error message 
# 

nohup ~/bin/EIG4.2/bin/smartpca -p m8.par > m8.Sout &

#######################################
# 5.3 to select individuals
cd ..

grep -v eigvals pca/m8.evec | \
 awk '{ if ( $2> 0.0044 && $2<0.0066 \
          && $3>-0.0044 && $3<0.0017 \
          && $4>-0.0054 && $4<0.0055 \
            ) print $1}' > t.ls

# use pc1,2,3 for the selection
# use the 5th and 6th 8 tiles for pc1
# use the 2nd and 3rd quartiles for pc2 and pc3

# 959 individuals
# 362 cases, 597 controls

sed -i 's/\:/ /' t.ls  
cp t.ls pca_ind.ls

# family ids and individual ids


#######################################
# 5.4 to generate the bed set

awk '{print $5}' bed/*_1.bim  | grep 0 | wc
# 6 SNPs with 0 (missing) for allele

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    awk '{print $2}' bed/${cohort}_1.bim | sort | uniq -d | wc
done
# no repeat of SNP ids in each cohort

# to get the common SNP list
rm t.ls

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    awk '{print $2}' bed/${cohort}_1.bim >> t.ls
done
# t.ls is now all snps in all cohorts

sort t.ls | uniq -c > t.out
# uniq -c will prefix lines by the number of occurrences
# takes about 10 mins

awk '{if ($1==8) print $2}' t.out > pca_snp.ls
# 2941089 SNPs 

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    plink --noweb --bfile bed/${cohort}_1 \
          --keep pca_ind.ls --extract pca_snp.ls \
          --make-bed --out ${cohort}_s
done

# to remove the SNPs with different locations in different cohorts
rm t.ls
for cohort in D1 D2 IR M3 NIT SLA UK; do 
    awk '{print $1,$2,$3,$4}' COR_s.bim > t1.map
    awk '{print $1,$2,$3,$4}' ${cohort}_s.bim > t2.map
    diff t1.map t2.map | grep ">" | awk '{print $3}' >> t.ls
done
sort t.ls | uniq > diff_loc_snp.ls
# 10484 SNPs 

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    plink --noweb --bfile ${cohort}_s \
          --exclude diff_loc_snp.ls \
          --make-bed --out ${cohort}_s
done
#2930605 SNPs left

# SNPs on different strands
rm t.ls 
for cohort in D1 D2 IR M3 NIT SLA UK; do 
    awk '{if ($6>$5){t=$6; $6=$5;$5=t} print $1,$2,$3,$4,$5,$6}' \
        COR_s.bim > t1.bim
    awk '{if ($6>$5){t=$6; $6=$5;$5=t} print $1,$2,$3,$4,$5,$6}' \
        ${cohort}_s.bim > t2.bim
    diff t1.bim t2.bim | grep ">" | awk '{print $3}' >> t.ls
done
sort t.ls | uniq > diff_strand_snp.ls
# 3310 SNPs 
for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    plink --noweb --bfile ${cohort}_s \
          --exclude diff_strand_snp.ls \
          --make-bed --out ${cohort}_s
done

# put bed bim fam filenames of other 7 cohorts into t.ls 
plink --noweb --bfile COR_s --merge-list t.ls --make-bed --out m8_s
# no warnings or error mesages...

# end up with 959 individuals and 2927295 SNPs 
# Total genotyping rate is 0.990583
awk '{if ($6==1) $6=2; if ($6==-9) $6=1; print $0}' m8_s.fam > t.out
mv t.out m8_s.fam 
# 362 cases and 597 controls 
mv m8_s.fam m8_s.bed m8_s.bim bed

#######################################
# 5.5 to get the grm file

~/bin/gcta64 --bfile bed/m8_s --autosome --make-grm --out grm/m8_s
# used about 12G memory 
# took 1.5 hours

awk '{print $1,$2,$6}' bed/m8_s.fam > m8_s.phen 
~/bin/gcta64 --grm grm/m8_s --pheno m8_s.phen --reml --prevalence 0.00005 
# var 0.112 SE 0.0274 
# size 959


################################################################################
# 6. the forest plot                                                           #
################################################################################

# input: 
#                       Var             SE         SIZE
# COR_1.hsq:V(G)/Vp_L   0.115020        0.043770    507 
# D1_1.hsq:V(G)/Vp_L    0.115042        0.027139    852 
# D2_1.hsq:V(G)/Vp_L    0.115019        0.008632   2599 
# IR_1.hsq:V(G)/Vp_L    0.115043        0.057633    420 
# M3_1.hsq:V(G)/Vp_L    0.118695        0.009268   2449 
# NIT_1.hsq:V(G)/Vp_L   0.126614        0.041567    653 
# SLA_1.hsq:V(G)/Vp_L   0.115205        0.007048   3536 
# UK_1.hsq:V(G)/Vp_L    0.141458        0.013874   2209 

#######################################
#in R:
library(rmeta)

# var 
v <- c(0.115020,0.115042,0.115019,0.115043,0.118695,0.126614,0.115205,0.141458)
# se 
se <- c(0.043770,0.027139,0.008632,0.057633,0.009268,0.041567,0.007048,0.013874)
# cohort sizes for weighting?
si <- c(   507,   852,  2599,   420,  2449,   653,  3536,  2209)

meta.summaries(v,se) 
# Fixed-effects meta-analysis
# Summary effect=0.119   95% CI (0.11, 0.127)
# meta.summaries(v,se,weights=si)
# with weighting, results didn't change much 

# Random-effects meta-analysis gave the same results

tabletext <- cbind(c('Cor','D1','D2','IR','M3','NIT','SLA','UK',NA,NA,'Summary',NA,NA))

# mean = Var and summary mean 
m <- c(0.1150,0.1150,0.1150,0.1150,0.1187,0.1266,0.1152,0.1415,NA,NA,0.119,NA,NA)
# lower = Var - 1.96*SE and summary lower 
l <- c(0.0292,0.0618,0.0981,0.0021,0.1005,0.0451,0.1014,0.1143,NA,NA,0.11,NA,NA)
# upper = Var + 1.96*SE and summary upper
u <- c(0.2008,0.1682,0.1319,0.2280,0.1369,0.2081,0.1290,0.1687,NA,NA,0.127,NA,NA)

pdf('meta_var.pdf')
forestplot(tabletext,m,l,u,is.summary=c(rep(FALSE,10),TRUE,FALSE,FALSE),,col=meta.colors(line="darkblue", summary="royalblue"),xticks=c(-0.1,0,0.1,0.2,0.3))
dev.off()

################################################################################
# 7. to investigate the prevalences, (10,7,5,3) * 10^-5                        #
################################################################################

# pv=0.00100
# pv=0.00050
# pv=0.00010
# pv=0.00007
# pv=0.00005
# pv=0.00003

pv=0.00003

# 4 PCs
for cohort in SLA D2; do 
    for p_thres in 0.5 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence ${pv} \
           --out ${cohort}_${p_thres}
    done
done

# 2 PCs
for cohort in M3 ; do
    for p_thres in 0.5 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence ${pv} \
           --out ${cohort}_${p_thres}
    done
done

# 1 PC
for cohort in IR NIT UK ; do
    for p_thres in 0.5 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} --qcovar ${cohort}.pc \
           --pheno ${cohort}_cc.phen --reml --prevalence ${pv} \
           --out ${cohort}_${p_thres}
    done
done

# no PC 
for cohort in D1 COR; do 
    for p_thres in 0.5 ; do
        ~/bin/gcta64 --grm grm/${cohort}_${p_thres} \
           --pheno ${cohort}_cc.phen --reml --prevalence ${pv} \
           --out ${cohort}_${p_thres}
    done
done

grep "V(G)/Vp_L"  *.hsq


