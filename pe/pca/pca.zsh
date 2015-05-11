################################################################################
#                         PCA of the PE GWAS data set                          #
################################################################################
 
# 1. Linkage disequilibrium based SNP pruning

# in the brc cluster 
# 22 jobs like
/share/bin/plink --noweb --bfile ~/scratch/dev/pca/pe10 \
  --chr 22 --hwe 1e-3 --maf 0.10 \
  --indep-pairwise 1500 150 0.2 --out pe10_22

# end with 71677 SNPs in pe10_s.ls

plink --noweb --bfile pe10 --extract pe10_s.ls --make-bed --out pe10_s

################################################################################
# 2. double check duplicates

# in the brc cluster
# 
split -l 1037 pe10_s.fam
cat xaa xab >  t0_1.fam
cat xaa xac >  t0_2.fam
...
cat xad xae >  t3_4.fam

plink --noweb --bfile pe10_s --keep t0_1.fam --make-bed --out p0_1&
plink --noweb --bfile pe10_s --keep t0_2.fam --make-bed --out p0_2&
...
plink --noweb --bfile pe10_s --keep t3_4.fam --make-bed --out p3_4&

~/bin/plink --noweb  --bfile ~/scratch/dev/pca/p0_1 \
 --genome  --min 0.05 \
 --out /home/klinbrc/scratch/dev/pca/p0_1
...
~/bin/plink --noweb  --bfile ~/scratch/dev/pca/p3_4 \
 --genome  --min 0.05 \
 --out /home/klinbrc/scratch/dev/pca/p3_4

# no duplications found 
# 1116 parent/offspring relationships, 36 are not presented in pe10_s.fam
# 1138 parent/offspring relationships in pe10_s.fam, 1080 (95%) are confirmed

################################################################################
# 3. merge with the Hapmap set

# I need GBR TSI IBS CHB JPT YRI 
# at least I need CEU TSI CHB YRI
# from ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/
# get the CEU, TSI, CHB, YRI from relationships_w_pops_121708.txt to t.fam
# 550 individuals
# have CEU TSI CHB YRI in t.in

grep -w -f t.in relationships_w_pops_121708.txt > t.fam

awk '{print $2}' pe10_s.bim > t.ls
plink --noweb --file hapmap3_r2_b36_fwd.consensus.qc.poly \
  --extract t.ls --keep t.fam --make-bed --out hapmap3

#  61806 SNPs 504 individuals

awk '{print $2}' hapmap3.bim > t.ls
plink --noweb --bfile pe10_s --extract t.ls --make-bed --out pe10_ss

plink --noweb --bfile pe10_ss --bmerge hapmap3.bed hapmap3.fam hapmap3.bim \
    --make-bed --out merge
# plink says I need a smaller data set or a bigger computer

awk '{print $2,$5,$6}' pe10_ss.bim > t1.ls
awk '{print $2,$5,$6}' hapmap3.bim > t2.ls

diff t1.ls t2.ls | grep ">" | awk '{print $2}' > t.ls

plink --noweb --bfile hapmap3 --flip t.ls --make-bed --out hapmap3_s

awk '{print $2,$5,$6}' pe10_ss.bim > t1.ls
awk '{print $2,$5,$6}' hapmap3_s.bim > t2.ls

diff t1.ls t2.ls | grep ">" | awk '{print $2}' > t.ls

plink --noweb --bfile pe10_ss --exclude t.ls --make-bed --out pe10_ss
plink --noweb --bfile hapmap3_s --exclude t.ls --make-bed --out hapmap3_s

diff hapmap3_s.bim pe10_ss.bim
# no difference found
plink --noweb --bfile pe10_ss --bmerge hapmap3_s.bim hapmap3_s.bed \
 hapmap3_s.fam --make-bed --out merg
#still says I need more memory

plink --noweb --bfile pe10_ss --recode --out pe10_ss
plink --noweb --bfile hapmap3_s --recode --out hapmap3_s

cat pe10_ss.ped hapmap3_s.ped > merge.ped
cp pe10_ss.map merge.map

plink --noweb --file merge --make-bed --out merge

# 51467 SNPs, 5685 individuals

# EIG wants individuals to be case or control, no missing population
# treat thre relatives and the 270 hapmap individuals at controls
awk '{if($6==-9 || $6==0)$6=1; print $0}' merge.fam > t.fam; mv t.fam merge.fam

################################################################################
# 4. PCA of the merged set

R --vanilla --slave --args stem=merge nsnpldregress=5 < EIGENSOFTplus_v12.txt
# smartpca wants 1G memory for 51467 SNPs, 5685 individuals

grep -w CEU relationships_w_pops_121708.txt | awk '{print $2}' > ceu.ls
grep -w TSI relationships_w_pops_121708.txt | awk '{print $2}' > tsi.ls
grep -w CHB relationships_w_pops_121708.txt | awk '{print $2}' > chb.ls
grep -w YRI relationships_w_pops_121708.txt | awk '{print $2}' > yri.ls
awk '{print $2}' pe10_s.fam > pe.ls

awk '{if ($2> $1*-1.8+0.004) print $0}' pe.dat | wc
# 246 individuals to be removed
# so who the hell are those guys?

grep -w PE_EDIN t.in    | awk '{print $1}' | sort -u > edin.ls
grep -w PE_HEID t.in    | awk '{print $1}' | sort -u > heid.ls
grep -w PE_HOLL t.in    | awk '{print $1}' | sort -u > holl.ls
grep -w PE_LOND t.in    | awk '{print $1}' | sort -u > lond.ls
grep -w PE_MUNI t.in    | awk '{print $1}' | sort -u > muni.ls
grep -w PE_PAMP t.in    | awk '{print $1}' | sort -u > pamp.ls
grep -w Valdecilla t.in | awk '{print $1}' | sort -u > vald.ls
grep -w WAFSS t.in      | awk '{print $1}' | sort -u > ause.ls

# there are 10 guys with no center attatched.
# WTCCCT534986 WTCCCT535009 WTCCCT535012 WTCCCT535057 WTCCCT535242
# WTCCCT535279 WTCCCT535290 WTCCCT535291 WTCCCT535338 WTCCCT536131

t1.py
----------------------------------------
#! /usr/bin/env python3

import sys

ifile=open(sys.argv[1])
ids=[];
for line in ifile:
    ids.append(line.strip())
ifile.close()

ifile=open('merge.evec')
ifile.readline()
for line in ifile:
    col=line.split()
    if col[0].split(':')[1] not in ids:
        continue;
    print(col[1],col[2])
ifile.close()
----------------------------------------

t1.py edin.ls > edin.dat
t1.py heid.ls > heid.dat
t1.py holl.ls > holl.dat
t1.py lond.ls > lond.dat
t1.py muni.ls > muni.dat
t1.py pamp.ls > pamp.dat
t1.py vald.ls > vald.dat
t1.py ause.ls > ause.dat

grace ????.dat &
mv merge_center_pca.png res
mv mergepca.png res

# now we can remove the 246 individuals, 
# forget about the other populations and do PCA on our PE set
t1.py 
----------------------------------------
#! /usr/bin/env python3

ids=[]
ifile=open('merge.evec')
ifile.readline()
for i in range(5181):
    line=ifile.readline()
    col=line.split()
    id=col[0].split(':')[1]
    x=float(col[1])
    y=float(col[2])
    if x*-1.8 + 0.004 < y :
        ids.append(id)
ifile.close()

ifile=open('pe10_s.fam')
for line in ifile:
    if line.split()[1] in ids:
        print(line[:-1])
ifile.close()
----------------------------------------
t1.py > t.fam
plink --noweb --bfile pe10_s --remove t.fam --make-bed --out pe11_s

awk '{if($6==-9 || $6==0)$6=1; print $0}' pe11_s.fam > t.fam;mv t.fam pe11_s.fam
# 4935 individuals, 71677 SNPs

################################################################################
# 5. PCA of the PE set

R --vanilla --slave --args stem=pe11_s nsnpldregress=5 < EIGENSOFTplus_v12.txt

# pc2 gave a lot outliers
# remove anyone beyond +-0.03
tail +2 pe11_s.evec | awk '{if ($3>0.03 || $3<-0.03) print $1}' | \
 sed 's/:/ /' > t.fam

# remove the 30 outliers in PC2
plink --noweb --bfile pe11_s --remove t.fam --make-bed --out pe12_s
# 4905 individuals
R --vanilla --slave --args stem=pe12_s nsnpldregress=5 < EIGENSOFTplus_v12.txt

tail +2 pe12_s.evec | awk '{if ($3>0.1) print $1}' | sed 's/:/ /' > t.fam
# remove 6 outliers in PC2
plink --noweb --bfile pe12_s --remove t.fam --make-bed --out pe13_s
# 4899 individuals

R --vanilla --slave --args stem=pe13_s nsnpldregress=5 < EIGENSOFTplus_v12.txt

tail +2 pe13_s.evec | awk '{if ($3>0.07) print $1}' | sed 's/:/ /' > t.fam
# remove 20 outliers in PC2
plink --noweb --bfile pe13_s --remove t.fam --make-bed --out pe14_s
# 4879 individuals

R --vanilla --slave --args stem=pe14_s nsnpldregress=5 < EIGENSOFTplus_v12.txt

tail +2 pe14_s.evec | awk '{if ($3<-0.1) print $1}' | sed 's/:/ /' > t.fam
# remove 5 outliers in PC2
plink --noweb --bfile pe14_s --remove t.fam --make-bed --out pe15_s
# 4874 individuals
R --vanilla --slave --args stem=pe15_s nsnpldregress=5 < EIGENSOFTplus_v12.txt

tail +2 pe15_s.evec | awk '{if ($3<-0.09 && $2 < 0.0) print $1}' | \
 sed 's/:/ /' > t.fam
tail +2 pe15_s.evec | awk '{if ($3>0.06) print $1}' | sed 's/:/ /' >> t.fam

plink --noweb --bfile pe15_s --remove t.fam --make-bed --out pe16_s
R --vanilla --slave --args stem=pe16_s nsnpldregress=5 < EIGENSOFTplus_v12.txt

################################################################################
# 6 double check the non-european individuals

plink --noweb --bfile merge --keep t.fam --make-bed --out ms
# the last 504 are ceu, tsi, chb, and yri individuals
awk '{print $2}' ms.bim > t.ls
plink --noweb --bfile pe16_s --extract t.ls --make-bed --out pe16_ss
diff pe16_ss.bim ms.bim | grep ">" | awk '{print $3}' > t.ls
plink --noweb --bfile ms --exclude t.ls --make-bed --out ms
plink --noweb --bfile pe16_ss --exclude t.ls --make-bed --out pe16_ss
diff pe16_ss.bim ms.bim 
# nothing
plink --noweb --bfile ms --bmerge pe16_ss.bim pe16_ss.bed pe16_ss.fam \
--make-bed --out mss
# complained about out of memory again

plink --noweb --bfile ms --recode --out ms
plink --noweb --bfile pe16_ss --recode --out pe16_ss
diff ms.map pe16_ss.map
cp ms.map mss.map
cat pe16_ss.ped ms.ped > mss.ped

plink --noweb --file mss --make-bed --out mss
# 50588 SNPs
# 4861+504=5365 individuals
R --vanilla --slave --args stem=mss nsnpldregress=5 < EIGENSOFTplus_v12.txt
tail +2 | tail -504 mss.evec | awk '{print $2,$3}' > yri.dat
tail +2  mss.evec | tail -504 | awk '{print $2,$3}' > yri.dat
tail +2 mss.evec | head -4861 | awk '{print $2,$3}' > pe.dat
# mss_center_pca.png in res/
# seems all caucasian

# from 5181 in pe10_s to 4861 in pe16_s
# 320 individuals removed
# mainly the 246 individuals removed being non-caucasian 

################################################################################
# 7 double check the ourliers 
# from the first logistic regression runs, it seems only the first three 
# covariates matter. 
# The first two PCs have no outliers. The third PC has five:
tail +2 pe16_s.evec | awk '{if($4>0.1)print $1}'
t.fam:
# F4003 WTCCCT526137
# F4003 WTCCCT526144
# F4003 WTCCCT524880
# F4003 WTCCCT526156
# F4003 WTCCCT524851

plink --noweb --bfile pe16_s --remove t.fam --make-bed --out pe17_s
R --vanilla --slave --args stem=pe17_s nsnpldregress=5 < EIGENSOFTplus_v12.txt

# from 5181 in pe10_s to 4856 individuals in pe17_s
# 325 individuals removed
# mainly the 246 individuals removed being non-caucasian 
################################################################################
#                    end of PCA of the PE GWAS data set                        #
################################################################################
