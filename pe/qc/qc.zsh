################################################################################
#                            QC of PE GWAS data                                #
################################################################################

# 1. Mendelian inheritance errors 

########################################
# 1.1 getting the Mendelian inheritance errors via pedstats tests
cp ../fam/pe3.??? .

sed 's/-9/0/' ../fam/pe3_part.fam > pe3.fam
awk '{if ($3!=0 && $4!=0) print $2,$3,$4}' pe3.fam | sort -k 1,1 > offspring.ls
# 364 offspring, father and mother triplets to check using pedstats
# 890661 SNPs

# pedstats wants both parents to be presented, or none at all.
# pedstats want 0 as missing for the disease phenotype

# out of memory on my 6G machine
# had to split the set into two
head -n 445330 pe3.bim | awk '{print $2}' > t1.ls
tail -n 445331 pe3.bim | awk '{print $2}' > t2.ls

plink --noweb --bfile pe3 --extract t1.ls --recode --out pe3_1
plink --noweb --bfile pe3 --extract t2.ls --recode --out pe3_2

awk '{if ($6==-9)$6=0;print $0}' pe4_1.ped > t.out; mv t.out pe4_1.ped
awk '{if ($6==-9)$6=0;print $0}' pe4_2.ped > t.out; mv t.out pe4_2.ped

# save 'A pe' into text file t.head
awk '{print "M",$2}' pe4_1.map > t.out; cat t.head t.out > pe4_1.dat
awk '{print "M",$2}' pe4_2.map > t.out; cat t.head t.out > pe4_2.dat

pedstats -p pe4_1.ped -d pe4_1.dat > ped_1.out
pedstats -p pe4_2.ped -d pe4_2.dat > ped_2.out

cat ped_1.out ped_2.out | grep -w Fam | awk '{print $1,$6}' > pedstats.out
# 877409 Mendelian inheritance errors detected
awk '{print $1}' pedstats.out| sort -g | uniq | wc
# 323620 SNPs (36.3%)
awk '{print $2}' pedstats.out| sort -g | uniq | wc
# all 352 triplets

awk '{print $2}' pedstats.out| sort | uniq -c | sort -g -r -k 1,1 > t.out
# on average, a triplets should have 919.4 SNP errors.
# the number is between 26999 and 684
awk '{print $1}' pedstats.out| sort | uniq -c | sort -g -r -k 1,1 > t.out
# for SNP, the number is between 168 and 0



########################################
# 1.2 missingness correlated to the number of Mendelian inheritance errors 

plink --noweb --bfile pe4 --missing --out pe4_missing

awk '{print $2}' pedstats.out | sort | uniq -c | awk '{print $2,$1}' | \
 sort -k 1,1 > t1.in 
# got the number of errors and the offspring ids

awk '{print $1}' offspring.ls > of.ls
grep -w -f of.ls pe4_missing.imiss | awk '{print $2,$6}' | sort -k 1,1 | \
 awk '{print $2}' > of.mis
paste t1.in of.mis > t.out
cp t.out t1.in 
# added the missingness of the offspring

rm of.ls of.mis

awk '{print $2,$1}' offspring.ls | sort -k 1,1 > fa.table
awk '{print $2}' offspring.ls | sort | uniq > fa.ls
grep -w -f fa.ls pe4_missing.imiss | awk '{print $2,$6}' | sort -k 1,1 > fa.mis 
join fa.table fa.mis | sort -k 2,2 | awk '{print $3}' >t.out; mv t.out fa.mis
# sort father missingness using the offspring id
paste t1.in fa.mis > t.out
cp t.out t1.in
# added the missingness of the father

rm fa.table fa.ls fa.mis

awk '{print $3,$1}' offspring.ls | sort -k 1,1 > mo.table
awk '{print $3}'  offspring.ls | sort | uniq > mo.ls
grep -w -f mo.ls pe4_missing.imiss | awk '{print $2,$6}' | sort -k 1,1 > mo.mis 
join mo.table mo.mis | sort -k 2,2 | awk '{print $3}' >t.out; mv t.out mo.mis
# sort mother missingness using the offspring id
paste t1.in mo.mis > t.out
cp t.out t1.in

rm mo.table mo.ls mo.mis

awk '{print $2,$3}' t1.in > t1.dat
awk '{print $2,$4}' t1.in > t2.dat
awk '{print $2,$5}' t1.in > t3.dat

t1.py
----------------------------------------
#! /usr/bin/env python3

ifile=open('t1.in')
for line in ifile:
    col=line.split()
    mis=max(float(col[2]),float(col[3]),float(col[4]))
    num=float(col[1])
    print(num,mis)
ifile.close()
----------------------------------------
t1.py > t.dat
grace t1.dat t.dat

# print to miss_vs_nme.png
# the number of Mendelian inheritance errors is related to the missingness

########################################
# 1.3 heterogeneous related to the number of Mendelian inheritance errors 

plink --noweb --bfile pe4 --het --out pe4_het

awk '{print $2}' pedstats.out | sort | uniq -c | awk '{print $2,$1}' | \
 sort -k 1,1 > t1.in 
# got the number of errors and the offspring ids

awk '{print $1}' offspring.ls > of.ls
grep -w -f of.ls pe4_het.het | awk '{print $2,$6}' | sort -k 1,1 | \
 awk '{print $2}' > of.het
paste t1.in of.het > t.out
cp t.out t1.in 
# added the heterogeneoty of the offspring

rm of.ls of.het

awk '{print $2,$3}' t1.in > t1.dat
awk '{print $2,$4}' t1.in > t2.dat
awk '{print $2,$5}' t1.in > t3.dat

# the comparison was saved into  het_vs_nme.png

########################################
# 1.4 remove SNPs and individuals

# to remove SNPs,
# use triplets with less than 5000 errors only,
awk '{if ($2> 5000) print $1}' t1.in > trpOff.ls
# 34 triplets needed to be removed

grep -w -v -f trpOff.ls pedstats.out | awk '{print $1}' | wc
# 482889 errors

grep -w -v -f trpOff.ls pedstats.out | awk '{print $1}' | sort | uniq | wc 
# instead of 323620 SNPs (36.3%), we now have 206542 (23.2%) SNPs to deal with 

grep -w -v -f trpOff.ls pedstats.out | awk '{print $1}' | sort | uniq -c | \
 sort -g -r -k 1,1 > t.in
 
# SNPs can fail in at most 152 triplets 

# the number of Mendelian inheritance errors is positively related to the 
# missingness of SNPs
# the correlation coefficient is 0.136201
# the t values are 62 and 290 for slope and constant 
# highly significient

awk '{if ($1>5)print $2}' t.in | awk '{if (rand()> 0.9932) print $0}' > t.ls
# select 100 from 13922
# 7 of the 100 are correctly genotyped
# 93 are not. 

awk '{if ($1==5)print $2}' t.in | awk '{if (rand()> 0.975585) print $0}' > t.ls
# 100 out of 4356 SNPs
# 19 of the 100 are correctly genotyped
# 81 are not

awk '{if ($1==4)print $2}' t.in | awk '{if (rand()> 0.98766) print $0}' > t.ls
# 100 out of 7882
# 20 are correct
# 80 are wrong

awk '{if ($1==1)print $2}' t.in | awk '{if (rand()> 0.99915) print $0}' > t.ls
# 100 out of 123366
# 27 correct
# 73 wrong

awk '{print $2}' ../fam/pe4.bim |  awk '{if (rand()> 0.999883) print $0}' > t.ls
# 100 out of 890661
# 51 correct
# 49 wrong
# should have done this ages ago...

########################################
# poisson distribution

t1.py
----------------------------------------
#! /usr/bin/env python3
from math import *

n=890661.0;
k=0;
la=482889/890661

k=0;p0=(la**k)*exp(-la)/factorial(k)
k=1;p1=(la**k)*exp(-la)/factorial(k)
k=2;p2=(la**k)*exp(-la)/factorial(k)
k=3;p3=(la**k)*exp(-la)/factorial(k)
k=4;p4=(la**k)*exp(-la)/factorial(k)
----------------------------------------
# p0+p1          = 0.89675
# p0+p1+p2       = 0.98221
# p0+p1+p2+p3    = 0.99766
# p0+p1+p2+p3+p4 = 0.99975

# it is 1% significent to cut all SNPs failed 4 or more times
# there should be 2087 if the background poisson is correct
# there are 26160, some 12 times more

awk '{if ($1>3) print $2}' t.in > t.ls
plink --noweb --bfile ../fam/pe4 --exclude t.ls --make-bed --out pe5

#5471 individuals 864501 SNPs

################################################################################
# 2. SNP missingness
plink --noweb --bfile pe5 --missing --out pe5_missing

t1.R
----------------------------------------
IMISS=read.table("pe5_missing.imiss", header=T, as.is=T)
LMISS=read.table("pe5_missing.lmiss", header=T, as.is=T)

png('pe5_missing.png')
oldpar=par(mfrow=c(1,2))
plot( (1:dim(IMISS)[1])/(dim(IMISS)[1]-1), sort(1-IMISS$F_MISS), 
 main="Ordered individual call rate", xlab="Quantile", ylab="Call Rate"); grid()
plot( (1:dim(LMISS)[1])/(dim(LMISS)[1]-1), sort(1-LMISS$F_MISS), 
 main="Ordered SNP coverage", xlab="Quantile", ylab="Coverage" ); grid()
par(oldpar)
dev.off()
----------------------------------------

# according to pe5_missing.png, use a threshold of 95% for SNPs

plink --noweb --bfile pe5 --geno 0.05 --make-bed --out pe6

# from 864501 SNPs to 852688 SNPs, removed 11813 SNPs

plink --noweb --bfile pe6 --test-missing --out pe6_phenomiss
plink --noweb --bfile pe6 --test-mishap  --out pe6_genomiss

plink --noweb --bfile pe6 --missing --out pe6_miss
plink --noweb --bfile pe6 --freq    --out pe6_freq

awk '{if (rand()>0.991) print $0}' pe6_genomiss.missing.hap | head


for ifile in  pe6_phenomiss.missing pe6_miss.lmiss pe6_freq.frq
do
    head -n 1 $ifile > t.head
    tail +2 $ifile | awk '{if (rand()>0.99) print $0}' > t.out
    cat t.head t.out > $ifile.s
done

ifile=pe6_genomiss.missing.hap
head -n 1 $ifile > t.head
tail +2 pe6_miss.lmiss.s | awk '{print $2}' > t.ls
grep -w -f t.ls $ifile > t.out  
# use a simple python script is way faster than using grep here.
cat t.head t.out > $ifile.s

t1.R
----------------------------------------
PMISS=read.table("pe6_phenomiss.missing.s", header=T, as.is=T)
GMISS=read.table("pe6_genomiss.missing.hap.s", header=T, as.is=T)
LMISS=read.table("pe6_miss.lmiss.s", header=T, as.is=T, row.names=2)
FREQ= read.table("pe6_freq.frq.s", header=T, as.is=T, row.names=2)

png('pe6_missing.png')

oldpar=par(mfrow=c(2,2))
OK = (LMISS$F_MISS>=0)&(FREQ$MAF>=0); OK[is.na(OK)]=FALSE
plot( LMISS$F_MISS[OK], PMISS$P[OK], main=NULL, 
 xlab="Missingness", ylab="P (miss vs pheno)" );
abline(lm(PMISS$P[OK]~LMISS$F_MISS[OK]),col='RED'); grid()

plot( FREQ$MAF[OK], PMISS$P[OK], main=NULL, 
 xlab="MAF", ylab="P (miss vs pheno)" ); 
abline(lm(PMISS$P[OK]~LMISS$F_MISS[OK]),col='RED'); grid()

OK=(GMISS$P>=0); OK[is.na(OK)]=FALSE
minP = tapply( GMISS$P[OK], as.factor(GMISS$SNP[OK]), min, na.rm=T )
minP=minP[is.finite(minP)]
OK=names(minP)
plot( LMISS[OK,"F_MISS"], minP, main=NULL, 
 xlab="Missingness", ylab="P (hap miss vs pheno)" ); 
abline(lm(minP~LMISS[OK,"F_MISS"]),col='RED'); grid()

plot( FREQ[OK,"MAF"], minP, main=NULL, 
 xlab="MAF", ylab="P (hap miss vs pheno)" );
abline(lm(minP~FREQ[OK,"MAF"]),col='RED'); grid()
par(oldpar)

dev.off()
----------------------------------------

################################################################################
# 3. SNP MAF
tail +2 pe6_freq.frq | awk '{print $5}' > t.out
listdis t.out > t.dat
grace t.dat

# print to pe6_maf.png
# use a 0.02 threshold

plink --noweb --bfile pe6 --maf 0.02 --make-bed --out pe7
# from 852688 to 707505, removed 145183 SNPs
# 5471 individuals

################################################################################
# 4. SNP Hardy-Weinberg Equilibrium
plink --noweb --bfile pe7 --hardy --out pe7_hwe
grep "O(HET)\|UNAFF" pe7_hwe.hwe > pe7_hwe_ctrl.hwe

t1.R
----------------------------------------
HWE=read.table("pe7_hwe_ctrl.hwe", header=T, as.is=T)
P = HWE$P[is.finite(HWE$P)]
n = length(P)
png('pe7_hwe_ctrl.png')
plot( qchisq((1:n)/(n+1),2), sort(-2*log(P)), 
 main="Q-Q plot of log(control HWE P-values)",
 xlab="Expected quantile", ylab="Observed quantile" )
lines( c(0,50), c(0,50) ); grid()
dev.off()
----------------------------------------
# this time the figure is slightly better, because we did better SNP cleaning
# use the same 1e-6 threshold

plink --noweb --bfile pe7 --hwe 1e-6 --make-bed --out pe8

# 1602 markers failed HWE test in cases
# 2409 markers failed HWE test in controls
# 705096 SNPs left

plink --noweb --bfile pe8 --hardy --out pe8_hwe
grep "O(HET)\|UNAFF" pe8_hwe.hwe > pe8_hwe_ctrl.hwe

# get the pe8_hwe_ctrl.png, it is still quite inflated.

################################################################################
# 5. inbreeding coefficients

# it is suggested in the plink websit to use --indep-pairwise and --indep SNPs
# use the last stratification set
awk '{print $2}' ../../PE/strat/pe8_s.bim > t.ls
# 71645 SNPs
plink --noweb --bfile pe8 --extract t.ls --make-bed --out pe8_s
# 70039 SNPs
plink --noweb --bfile pe8_s --het --out pe8_het

t1.R
----------------------------------------
HET=read.table("pe8_het.het", header=T, as.is=T)
H = (HET$N.NM.-HET$O.HOM.)/HET$N.NM.
png('pe8_het.png')
hist(HET$F,50)
dev.off()
----------------------------------------

# set -0.07 <= F <= 0.07...
tail +2 pe8_het.het |awk '{if ($6>0.07 || $6 < -0.07) print $2}' > t.ls
# 82 individuals

t1.py
----------------------------------------
#! /usr/bin/env python3

ifile=open('t.ls')
ind=[];
for line in ifile:
    ind.append(line.strip())
ifile.close()

ifile=open('pe8.fam')
for line in ifile:
    col=line.split()
    if col[1] not in ind:
        continue;
    print(line[:-1])
ifile.close()
----------------------------------------
t1.py > t.fam

plink --noweb --bfile pe8 --remove t.fam --make-bed --out pe9
# from 5471 individuals to 5389
# 705096 SNPs

################################################################################
# 6. individual missingness

plink --noweb --bfile pe9 --missing --out pe9_missing

# got the pe9_missing.png

# if I applied the 0.02 threshold again, I end up with 5014 individuals,
# 2 more than 5012 we had last time before stratification

plink --noweb --bfile pe9 --mind 0.03 --make-bed --out pe10 
# from 5389 to 5181 individuals, removed 208
# 705096 SNPs

