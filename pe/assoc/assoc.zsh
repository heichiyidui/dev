################################################################################
#            association studies of the PE GWAS data set                       #
################################################################################

1. the SNPs manually removed in the last round

awk '{print $2}' ../../PE/strat/pe5.bim | sort > t1.ls
awk '{print $2}' ../../PE/assoc/plink/pe11.bim | sort > t2.ls
diff t1.ls t2.ls | grep "<" | awk '{print $2}' > to_double_check.ls

# Evoker to pdf
awk '{print $1 "\t-1"}' to_double_check.ls > to_double_check.ls.scores
# to pdf, pdftohtml then

t1.py 
----------------------------------------
#! /usr/bin/env python3

ifile=open('thtml/alls.html')
for line in ifile:
    pid='';
    if line.startswith('<A name='):
        pid=line.split('\"')[1]
        line=ifile.readline()
        print('mv '+pid+' '+line.split('<')[0]+'.png')
ifile.close()
----------------------------------------

# 4260 SNPs removed in the last round
# 3177 of them in pe10.bim in this round
# saved in to_double_check.ls
# found 752 SNPs to be reconsidered 
# 2425 SNPs to be deleted. Hopefully I'll never see them again. 

awk '{if ($2==-1) print $1}' res/to_double_check.ls.scores > t.ls
plink --noweb --bfile ../qc/pe10 --keep ../pca/pe16_s.fam \
 --exclude t.ls --maf 0.02 --make-bed --out pe17
# 4861 individuals
# 696631 SNPs
# 6040 SNPs removed for maf < 0.02
# 2425 SNPs removed for bad genotyping

################################################################################
# 2. number of PCs

vi t.head
FID IID COV1 COV2 COV3 COV4 COV5 COV6 COV7 COV8 COV9 COV10

tail +2 ../pca/pe16_s.evec | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' \
 | sed 's/:/ /' > t.out
cat t.head t.out > t.cov

########################################

plink --noweb --bfile pe17 --logistic --out cov0
plink --noweb --bfile pe17 --logistic --covar t.cov --covar-number 1 \
 --out cov1
plink --noweb --bfile pe17 --logistic --covar t.cov --covar-number 1-2 \
 --out cov2
#...
plink --noweb --bfile pe17 --logistic --covar t.cov --covar-number 1-9 \
 --out cov9
 
# from http://gettinggeneticsdone.blogspot.com/2011/04/ \
# annotated-manhattan-plots-and-qq-plots.html

awk '{if ($5 == "ADD" || $5 =="TEST") print $1,$2,$3,$9}' cov4.assoc.logistic \
 > cov4.dat
t1.R
----------------------------------------
source("http://www.StephenTurner.us/qqman.r")

results <- read.table("cov0.dat",T);
png(filename='cov0_man.png',width=1600,height=600);manhattan(results);
dev.off();
png('cov0_qq.png',width=1200,height=1200);qq(results$P);dev.off();
----------------------------------------

# it seems 3 covariates should be enough

grep -w -v CHR cov?.dat | awk '{if ($4<0.001) print $2}' | \
 sort -u > cov_snp_chk.ls
awk '{print $1 "\t-1"}' cov_snp_chk.ls > cov_snp_chk.ls.scores
# Evoker to pdf, pdftohtml, t1.py etc. 
# out of the 4342 SNPs, 1037 are to be removed.


################################################################################
# 3. Go back to PCA, correct the first three PCs.
# removed the 5 individuals of the F4003 family

plink --noweb --bfile pe17 --remove ../pca/t.fam --make-bed --out pe17

tail +2 ../pca/pe17_s.evec | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' \
 | sed 's/:/ /' > t.out
cat t.head t.out > t.cov

# remove the 1037 badly genotyped SNPs found in the first round of cov runs

awk '{if ($2==-1)print $1}' cov_snp_chk.ls.scores > t.ls
plink --noweb --bfile pe17 --exclude t.ls --make-bed --out pe17
# 4856 individuals
# from 696631 SNPs to 695594 SNPs, 1037 removed

tail +2 ../pca/pe17_s.evec | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' \
 | sed 's/:/ /' > t.out
cat t.head t.out > t.cov

################################################################################
# 4. second round of COV runs

# 3 PCs are needed

grep -w -v CHR cov?.dat | awk '{if ($4<0.001) print $2}' | \
 sort -u > cov_snp_chk2.ls
# 3352 SNPs
awk '{print $1 "\t-1"}' cov_snp_chk2.ls > cov_snp_chk2.ls.scores

# 3352 SNPs to be manually checked.
# 59 SNPs removed
awk '{if ($2==-1)print $1}' cov_snp_chk2.ls.scores > t.ls
plink --noweb --bfile pe17 --exclude t.ls --make-bed --out pe17 

##########################
# we have 2 manual checks of the top COV SNPs. For the SNPs past both checks,
# I don't want to see them any more.
awk '{if ($2==1) print $1}' res/cov*.scores | sort | uniq -c | \
 awk '{if ($1==2) print $2}' > res/double_pass.ls
# 3257 SNPs

################################################################################
# 5. unphased run

# in the unphased directory
cp ../pe17.??? .
awk '{if ($6==-9)$6=1; print $0}' pe17.fam > t.out
cp t.out pe17.fam
# now the relatives can be regarded as controls

# use only 3 covariates
grep -w -v FID ../t.cov | awk '{print $3,$4,$5}' > t.cov

plink --noweb --bfile pe17 --recode --out pe17

awk '{print $2}' pe17.bim > t.ls
# 695535 SNPs
split -l 6956 t.ls
# mv the names from xaa ... xdv to t00 ... t99

----------------------------------------
#!~/bin/sh
#$-S ~/bin/sh
for i in ~/scratch/dev/assoc/unphased/t??
do
    ~/bin/plink --noweb --bfile ~/scratch/dev/assoc/unphased/pe17 \
    --extract $i --recode --out $i
done
----------------------------------------

# paste t00.ped t.cov > t.ped; mv t.ped t00.3cov.ped

for i in t??
do
    paste $i.ped t.cov > t.ped; mv t.ped $i.3cov.ped
done

# creat two files: t.head and t.end
# t.head: "A      case"
# t.end:
# C       cov1
# C       cov2
# C       cov3

# awk '{print "M\t" $2}' t00.map > t.map; 
# cat t.head t.map t.end > t00.3cov.data
for i in t??
do
    awk '{print "M\t" $2}' $i.map > t.map
    cat t.head t.map t.end > $i.3cov.data
done

# prepared the input files

cf3_00.sh
----------------------------------------
#!/bin/sh
#$ -S /bin/sh
#$ -o  /home/klinbrc/scratch/dev/assoc/unphased
#$ -e  /home/klinbrc/scratch/dev/assoc/unphased

~/bin/unphased ~/scratch/dev/assoc/unphased/t00.3cov.ped \
         -data ~/scratch/dev/assoc/unphased/t00.3cov.data \
         -confounder cov1 cov2 cov3  > cf3_pe17_00.out
----------------------------------------

for i in cf3_??.sh
do 
    qsub -q long.q,short.q $i
done

cat cf3_*.out > up.out
grep Likeli up.out | awk '{print $NF}' > p.ls
awk '{print $1,$2,$4}' pe17.bim > t.in
paste t.in p.ls > t.out
# t.head : CHR SNP BP P
cat t.head t.out > up.table

t1.R:
----------------------------------------
source("http://www.StephenTurner.us/qqman.r")

results <- read.table("up.table",T);
png(filename='up_3_man.png',width=1600,height=600);manhattan(results);dev.off();
png('up_3_qq.png',width=1200,height=1200);qq(results$P);dev.off();
----------------------------------------
grep -w -v CHR up.table | awk '{if ($4<0.001) print $2}' | \
 sort -u > up_snp_chk.ls
# 765 SNPs to check
grep -w -v -f res/double_pass.ls up_snp_chk.ls > t.ls
mv t.ls up_snp_chk.ls
# 412 to check, 353 past twice already
awk '{print $1 "\t-1"}' up_snp_chk.ls > up_snp_chk.ls.scores
# after Evoker check, 341 left, 71 deleted.

awk '{if ($2==-1) print $1}' res/up_snp_chk.ls.scores > t.ls
grep -w -v -f t.ls up.table > t.out
mv t.out up.table 

plink --bfile pe17 --noweb --exclude t.ls --make-bed --out pe17
# 695464 SNPs left
# 4856 individuals

# no significent hits...
################################################################################
