# the phenotype file: data/Quant_traits_resid.csv
# there are four phenotypes, WBV, LVV, P3A, P3L
# Use the residues. Ranking is no much different.

################################################################################
# 1. prepare plink binary files                                                #
################################################################################

awk '{print $1,$2}' ../replica/pe19.fam > pe19_12.fam

################################################################################
# WBV
awk -F "," '{if ($17!="NA") print $1}' data/Quant_traits_resid.csv | \
 grep -v SANGERID > t.ls
# 1345 sangerids

grep -f t.ls pe19_12.fam > t.fam
# 732 individuals

plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out pewbv

#t1.py:
#---------------------------------------
#! /usr/bin/env python3
import sys

col=int(sys.argv[2])

ifile=open(sys.argv[1])
pheno={}
ifile.readline()
for line in ifile:
    if line.split(',')[col] != 'NA':
        pheno[line.split(',')[0]]=line.split(',')[col]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()
#---------------------------------------
t1.py data/Quant_traits_resid.csv 16 > t.phe

paste pewbv.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pewbv.fam

################################################################################
# LVV
awk -F "," '{if ($21!="NA") print $1}' data/Quant_traits_resid.csv | \
 grep -v SANGERID > t.ls
# 1345 sangerids

grep -f t.ls pe19_12.fam > t.fam
# 731 individuals

plink --noweb --bfile ../replica/pe19  --keep t.fam --make-bed --out pelvv

t1.py data/Quant_traits_resid.csv 20 > t.phe

paste pelvv.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pelvv.fam

################################################################################
# p3a
awk -F "," '{if ($25!="NA") print $1}' data/Quant_traits_resid.csv | \
 grep -v SANGERID > t.ls
# 1041 sanger ids
grep -f t.ls pe19_12.fam > t.fam
# 497 individuals
plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out pep3a

t1.py data/Quant_traits_resid.csv 24 > t.phe

paste pep3a.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pep3a.fam

################################################################################
# p3l
awk -F "," '{if ($29!="NA") print $1}' data/Quant_traits_resid.csv | \
 grep -v SANGERID > t.ls
# 1053 sanger ids
grep -f t.ls pe19_12.fam > t.fam
# 506 individuals
plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out pep3l

t1.py data/Quant_traits_resid.csv 28 > t.phe

paste pep3l.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pep3l.fam

################################################################################
# the T1, IMM, DEL phenotypes
awk -F","  '{if ($2!="")print $1}' \
 data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv | \
 grep -v SANGERID > t.ls
 
# 3513 individuals
grep -f t.ls pe19_12.fam > t.fam
# 2391 individuals

plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out pet10

t1.py data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv 1 > t.phe

paste pet10.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pet10.fam

#######################################
awk -F","  '{if ($3!="")print $1}' \
 data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv | \
 grep -v SANGERID > t.ls
# 3532 individuals

grep -f t.ls pe19_12.fam > t.fam
# 2400 individuals

plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out peimm

t1.py data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv 2 > t.phe

paste peimm.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out peimm.fam

#######################################
awk -F","  '{if ($4!="")print $1}' \
 data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv | \
 grep -v SANGERID > t.ls
# 3498 individuals

grep -f t.ls pe19_12.fam > t.fam
# 2378 individuals

plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out pedel

t1.py data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv 3 > t.phe

paste pedel.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pedel.fam

#######################################
# block pheno
awk -F","  '{if ($2!="")print $1}' data/BLOCK-DIGIT_23FEB12_4Kuang.csv \
 | grep -v SANGERID > t.ls
# 5516 individuals
grep -f t.ls pe19_12.fam > t.fam
# 3094 individuals

plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out peblo
t1.py data/BLOCK-DIGIT_23FEB12_4Kuang.csv 1 > t.phe

paste peblo.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out peblo.fam

#######################################
# digit pheno
awk -F","  '{if ($3!="")print $1}' data/BLOCK-DIGIT_23FEB12_4Kuang.csv \
 | grep -v SANGERID > t.ls
# 3131 individuals
grep -f t.ls pe19_12.fam > t.fam
# 1440 individuals

plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out pedig

t1.py data/BLOCK-DIGIT_23FEB12_4Kuang.csv 2 > t.phe
paste pedig.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pedig.fam

################################################################################
# 2. generate ped and map files                                                #
################################################################################
for i in blo del dig imm lvv p3a p3l t10 wbv
do 
    plink --noweb --bfile bed/pe$i --recode --out pe$i
     
    awk '{print $2}' bed/pe$i.fam > t.ls
    grep -f t.ls pe17_s.evec | awk '{print $2,$3,$4}' > t.cov

    paste pe$i.ped t.cov > t.ped; mv t.ped pe$i.ped

    awk '{print "M\t" $2}' pe$i.map > t.map
    cat t.head t.map t.end > pe$i.data
done

################################################################################
# 3. run UNPHASED                                                              #
################################################################################

#!/bin/sh
#$ -S /bin/sh
#$ -o  /scratch/home/klinbrc/dev/pheno
#$ -e  /scratch/home/klinbrc/dev/pheno

~/bin/unphased /scratch/home/klinbrc/dev/pheno/peblo.ped  \
         -data /scratch/home/klinbrc/dev/pheno/peblo.data  > peblo.out


for i in t?.sh
do
    qsub -q long.q,short.q -h_vmem 10G $i
done
    
# or with the three cov

#!/bin/sh
#$ -S /bin/sh
#$ -o  /scratch/home/klinbrc/dev/pheno
#$ -e  /scratch/home/klinbrc/dev/pheno

~/bin/unphased /scratch/home/klinbrc/dev/pheno/peblo.ped  \
         -data /scratch/home/klinbrc/dev/pheno/peblo.data \
         -confounder cov1 cov2 cov3  > peblo_cov3.out


################################################################################
# repeat that for the ranking-normaled pheno types                             #
################################################################################
# WBV
# awk -F "," '{if ($19!="NA") print $1}' data/Quant_traits_resid.csv | \
# grep -v SANGERID > t.ls
# 1353 sangerids

# grep -f t.ls pe19_12.fam > t.fam
# 736 individuals

# plink --noweb --bfile ../replica/pe19 --keep t.fam --make-bed --out pewbv
#t1.py data/Quant_traits_resid.csv 18 > t.phe
#paste pewbv.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
#mv t.out pewbv.fam
# ... the ranking results are about the same
################################################################################


################################################################################
# rank normalized pheno run finished                                           #
################################################################################

awk '{print $2}' pe18.bim > snp.ls

rm hits.ls
grep p-value res/rankn/pelvv.out | awk '{print $NF}' > t.out; 
paste snp.ls t.out > t.table; awk '{if ($2<0.0001) print $1}' t.table >> hits.ls
grep p-value res/rankn/pep3a.out | awk '{print $NF}' > t.out; 
paste snp.ls t.out > t.table; awk '{if ($2<0.0001) print $1}' t.table >> hits.ls
grep p-value res/rankn/pep3l.out | awk '{print $NF}' > t.out; 
paste snp.ls t.out > t.table; awk '{if ($2<0.0001) print $1}' t.table >> hits.ls
grep p-value res/rankn/pewbv.out | awk '{print $NF}' > t.out; 
paste snp.ls t.out > t.table; awk '{if ($2<0.0001) print $1}' t.table >> hits.ls
grep p-value res/resid/pelvv.out | awk '{print $NF}' > t.out; 
paste snp.ls t.out > t.table; awk '{if ($2<0.0001) print $1}' t.table >> hits.ls
grep p-value res/resid/pep3a.out | awk '{print $NF}' > t.out; 
paste snp.ls t.out > t.table; awk '{if ($2<0.0001) print $1}' t.table >> hits.ls
grep p-value res/resid/pep3l.out | awk '{print $NF}' > t.out; 
paste snp.ls t.out > t.table; awk '{if ($2<0.0001) print $1}' t.table >> hits.ls
grep p-value res/resid/pewbv.out | awk '{print $NF}' > t.out; 
paste snp.ls t.out > t.table; awk '{if ($2<0.0001) print $1}' t.table >> hits.ls

sort -u hits.ls > t.out ; mv t.out hits.ls

# 1932 hits
# Evoker to generate pics
# alls.html
grep "A name=" pics/alls.html| awk -F"\"" '{print $2}' > pic.ls
grep "rs" pics/alls.html | sed 's/<br\/>/.png/' > tg.ls
paste pic.ls tg.ls > t.in
awk '{print "mv " $0}' t.in > t.sh
cd pics 
source ../t.sh

ls pics/rs* | sed 's/pics\///' | sed 's/.png//' | sort > t.ls
# 1675 pictures left
sort hits.ls > t2.ls
# 1932 pictures 
diff t.ls t2.ls | grep ">" | awk '{print $2}' > to_remove.ls
# 257 snps to remove

################################################################################
# to get the manhattan plots                                                   #
################################################################################
t.head :
SNP CHR BP P

# snp, chr, bp
awk '{print $2,$1,$4}' pe18.bim > t.table

grep Likeli res/rankn/pelvv.out | awk '{print $NF}' > t.ls; 
paste t.table t.ls > t.out; 
cat t.head t.out | grep -w -v -f to_remove.ls > rankn_lvv.table

# do it for all 8 files 

R

source("http://www.StephenTurner.us/qqman.r")

results <- read.table("rankn_lvv.table",T);
png(filename='rankn_lvv_man.png',width=1600,height=600);
manhattan(results);dev.off();
png('rankn_lvv_qq.png',width=1200,height=1200);
qq(results$P);dev.off();

################################################################################
# found problems with QQ plots, they are all inflated. Run without cov         #
################################################################################

# phenotype RAVLT
awk  '{print $1}' data/RAVLT_TRIAL1_4Kuang_n_Cathryn_3Feb12final.csv | \
 grep -v SANGERID > t.ls
# 3072 individuals
grep -f t.ls pe18_12.fam > t.fam
# 2402 individuals

plink --noweb --bfile pe18 --keep t.fam --make-bed --out perav
t1.py
----------------------------------------
#! /usr/bin/env python3

ifile=open('data/RAVLT_TRIAL1_4Kuang_n_Cathryn_3Feb12final.csv')
pheno={}
ifile.readline()
for line in ifile:
    pheno[line.split()[0]]=line.split()[7]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()
----------------------------------------
t1.py > t.phe

paste perav.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out perav.fam

plink --noweb --bfile perav --recode --out perav

t.head: "T      rav"
awk '{print "M\t" $2}' perav.map > t.map
cat t.head t.map > perav.data
~/bin/unphased perav.ped -data perav.data > perav.out
################################################################################
# to get 10 subsets
awk '{print $2}' perav.bim > snp.ls
# 695464 snps
split -l 69547 snp.ls

plink --noweb --bfile perav --extract xaa --recode --out perav00
plink --noweb --bfile perav --extract xab --recode --out perav01
plink --noweb --bfile perav --extract xac --recode --out perav02
plink --noweb --bfile perav --extract xad --recode --out perav03
plink --noweb --bfile perav --extract xae --recode --out perav04
plink --noweb --bfile perav --extract xaf --recode --out perav05
plink --noweb --bfile perav --extract xag --recode --out perav06
plink --noweb --bfile perav --extract xah --recode --out perav07
plink --noweb --bfile perav --extract xai --recode --out perav08
plink --noweb --bfile perav --extract xaj --recode --out perav09

awk '{print "M\t" $2}' perav00.map > t.map; cat t.head t.map > perav00.data
awk '{print "M\t" $2}' perav01.map > t.map; cat t.head t.map > perav01.data
awk '{print "M\t" $2}' perav02.map > t.map; cat t.head t.map > perav02.data
awk '{print "M\t" $2}' perav03.map > t.map; cat t.head t.map > perav03.data
awk '{print "M\t" $2}' perav04.map > t.map; cat t.head t.map > perav04.data
awk '{print "M\t" $2}' perav05.map > t.map; cat t.head t.map > perav05.data
awk '{print "M\t" $2}' perav06.map > t.map; cat t.head t.map > perav06.data
awk '{print "M\t" $2}' perav07.map > t.map; cat t.head t.map > perav07.data
awk '{print "M\t" $2}' perav08.map > t.map; cat t.head t.map > perav08.data
awk '{print "M\t" $2}' perav09.map > t.map; cat t.head t.map > perav09.data

~/bin/unphased perav00.ped -data perav00.data > perav00.out
~/bin/unphased perav01.ped -data perav01.data > perav01.out
~/bin/unphased perav02.ped -data perav02.data > perav02.out
~/bin/unphased perav03.ped -data perav03.data > perav03.out
~/bin/unphased perav04.ped -data perav04.data > perav04.out
~/bin/unphased perav05.ped -data perav05.data > perav05.out
~/bin/unphased perav06.ped -data perav06.data > perav06.out
~/bin/unphased perav07.ped -data perav07.data > perav07.out
~/bin/unphased perav08.ped -data perav08.data > perav08.out
~/bin/unphased perav09.ped -data perav09.data > perav09.out

# cat them together, get perav.out

################################################################################
# get the R plots

################################################################################
# the T1, IMM, DEL phenotypes

awk -F","  '{if ($2!="")print $1}' \
 data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv | \
 grep -v SANGERID > t.ls
# 3513 individuals
grep -f t.ls pe18_12.fam > t.fam
# 2404 individuals

plink --noweb --bfile pe18 --keep t.fam --make-bed --out pet10
t1.py
----------------------------------------
#! /usr/bin/env python3

ifile=open('data/RAY ALL VARIABLES - 10thFeb-2Kuang-1.csv')
pheno={}
ifile.readline()
for line in ifile:
    pheno[line.split(',')[0]]=line.split(',')[1]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()

----------------------------------------
t1.py > t.phe

paste pet10.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pet10.fam

plink --noweb --bfile pet10 --recode --out pet10

t.head: "T      t1"
awk '{print "M\t" $2}' pet10.map > t.map
cat t.head t.map > pet10.data
~/bin/unphased pet10.ped -data pet10.data > pet10.out
################################################################################
awk -F","  '{if ($3!="")print $1}' \
 data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv | \
 grep -v SANGERID > t.ls
# 3532 individuals
grep -f t.ls pe18_12.fam > t.fam
# 2413 individuals

plink --noweb --bfile pe18 --keep t.fam --make-bed --out peimm
t1.py
----------------------------------------
#! /usr/bin/env python3

ifile=open('data/RAY ALL VARIABLES - 10thFeb-2Kuang-1.csv')
pheno={}
ifile.readline()
for line in ifile:
    pheno[line.split(',')[0]]=line.split(',')[2]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()

----------------------------------------
t1.py > t.phe

paste peimm.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out peimm.fam

plink --noweb --bfile peimm --recode --out peimm

t.head: "T      imm"
awk '{print "M\t" $2}' peimm.map > t.map
cat t.head t.map > peimm.data
~/bin/unphased peimm.ped -data peimm.data > peimm.out

################################################################################
awk -F","  '{if ($4!="")print $1}' \
 data/RAY\ ALL\ VARIABLES\ -\ 10thFeb-2Kuang-1.csv | \
 grep -v SANGERID > t.ls
# 3498 individuals
grep -f t.ls pe18_12.fam > t.fam
# 2391 individuals

plink --noweb --bfile pe18 --keep t.fam --make-bed --out pedel
t1.py
----------------------------------------
#! /usr/bin/env python3

ifile=open('data/RAY ALL VARIABLES - 10thFeb-2Kuang-1.csv')
pheno={}
ifile.readline()
for line in ifile:
    pheno[line.split(',')[0]]=line.split(',')[3]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()

----------------------------------------
t1.py > t.phe

paste pedel.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pedel.fam

plink --noweb --bfile pedel --recode --out pedel

t.head: "T      del"
awk '{print "M\t" $2}' pedel.map > t.map
cat t.head t.map > pedel.data
~/bin/unphased pedel.ped -data pedel.data > pedel.out
################################################################################

# to get 10 subsets
awk '{print $2}' pedel.bim > snp.ls
# 695464 snps
split -l 69547 snp.ls

plink --noweb --bfile pedel --extract xaa --recode --out pedel00
plink --noweb --bfile pedel --extract xab --recode --out pedel01
plink --noweb --bfile pedel --extract xac --recode --out pedel02
plink --noweb --bfile pedel --extract xad --recode --out pedel03
plink --noweb --bfile pedel --extract xae --recode --out pedel04
plink --noweb --bfile pedel --extract xaf --recode --out pedel05
plink --noweb --bfile pedel --extract xag --recode --out pedel06
plink --noweb --bfile pedel --extract xah --recode --out pedel07
plink --noweb --bfile pedel --extract xai --recode --out pedel08
plink --noweb --bfile pedel --extract xaj --recode --out pedel09

awk '{print "M\t" $2}' pedel00.map > t.map; cat t.head t.map > pedel00.data
awk '{print "M\t" $2}' pedel01.map > t.map; cat t.head t.map > pedel01.data
awk '{print "M\t" $2}' pedel02.map > t.map; cat t.head t.map > pedel02.data
awk '{print "M\t" $2}' pedel03.map > t.map; cat t.head t.map > pedel03.data
awk '{print "M\t" $2}' pedel04.map > t.map; cat t.head t.map > pedel04.data
awk '{print "M\t" $2}' pedel05.map > t.map; cat t.head t.map > pedel05.data
awk '{print "M\t" $2}' pedel06.map > t.map; cat t.head t.map > pedel06.data
awk '{print "M\t" $2}' pedel07.map > t.map; cat t.head t.map > pedel07.data
awk '{print "M\t" $2}' pedel08.map > t.map; cat t.head t.map > pedel08.data
awk '{print "M\t" $2}' pedel09.map > t.map; cat t.head t.map > pedel09.data

~/bin/unphased pedel00.ped -data pedel00.data > pedel00.out
~/bin/unphased pedel01.ped -data pedel01.data > pedel01.out
~/bin/unphased pedel02.ped -data pedel02.data > pedel02.out
~/bin/unphased pedel03.ped -data pedel03.data > pedel03.out
~/bin/unphased pedel04.ped -data pedel04.data > pedel04.out
~/bin/unphased pedel05.ped -data pedel05.data > pedel05.out
~/bin/unphased pedel06.ped -data pedel06.data > pedel06.out
~/bin/unphased pedel07.ped -data pedel07.data > pedel07.out
~/bin/unphased pedel08.ped -data pedel08.data > pedel08.out
~/bin/unphased pedel09.ped -data pedel09.data > pedel09.out
################################################################################

################################################################################
# block pheno
awk -F","  '{if ($2!="")print $1}' data/BLOCK-DIGIT_23FEB12_4Kuang.csv \
 | grep -v SANGERID > t.ls
# 5516 individuals
grep -f t.ls pe18_12.fam > t.fam
# 3105 individuals

plink --noweb --bfile pe18 --keep t.fam --make-bed --out peblo
t1.py
----------------------------------------
#! /usr/bin/env python3

ifile=open('data/BLOCK-DIGIT_23FEB12_4Kuang.csv')
pheno={}
ifile.readline()
for line in ifile:
    pheno[line.split(',')[0]]=line.split(',')[1]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()

----------------------------------------
t1.py > t.phe

paste peblo.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out peblo.fam

plink --noweb --bfile peblo --recode --out peblo

t.head: "T      blo"
awk '{print "M\t" $2}' peblo.map > t.map
cat t.head t.map > peblo.data
~/bin/unphased peblo.ped -data peblo.data > peblo.out
################################################################################

# to get 10 subsets
awk '{print $2}' peblo.bim > snp.ls
# 695464 snps
split -l 69547 snp.ls

plink --noweb --bfile peblo --extract xaa --recode --out peblo00
plink --noweb --bfile peblo --extract xab --recode --out peblo01
plink --noweb --bfile peblo --extract xac --recode --out peblo02
plink --noweb --bfile peblo --extract xad --recode --out peblo03
plink --noweb --bfile peblo --extract xae --recode --out peblo04
plink --noweb --bfile peblo --extract xaf --recode --out peblo05
plink --noweb --bfile peblo --extract xag --recode --out peblo06
plink --noweb --bfile peblo --extract xah --recode --out peblo07
plink --noweb --bfile peblo --extract xai --recode --out peblo08
plink --noweb --bfile peblo --extract xaj --recode --out peblo09

awk '{print "M\t" $2}' peblo00.map > t.map; cat t.head t.map > peblo00.data
awk '{print "M\t" $2}' peblo01.map > t.map; cat t.head t.map > peblo01.data
awk '{print "M\t" $2}' peblo02.map > t.map; cat t.head t.map > peblo02.data
awk '{print "M\t" $2}' peblo03.map > t.map; cat t.head t.map > peblo03.data
awk '{print "M\t" $2}' peblo04.map > t.map; cat t.head t.map > peblo04.data
awk '{print "M\t" $2}' peblo05.map > t.map; cat t.head t.map > peblo05.data
awk '{print "M\t" $2}' peblo06.map > t.map; cat t.head t.map > peblo06.data
awk '{print "M\t" $2}' peblo07.map > t.map; cat t.head t.map > peblo07.data
awk '{print "M\t" $2}' peblo08.map > t.map; cat t.head t.map > peblo08.data
awk '{print "M\t" $2}' peblo09.map > t.map; cat t.head t.map > peblo09.data

~/bin/unphased peblo00.ped -data peblo00.data > peblo00.out
~/bin/unphased peblo01.ped -data peblo01.data > peblo01.out
~/bin/unphased peblo02.ped -data peblo02.data > peblo02.out
~/bin/unphased peblo03.ped -data peblo03.data > peblo03.out
~/bin/unphased peblo04.ped -data peblo04.data > peblo04.out
~/bin/unphased peblo05.ped -data peblo05.data > peblo05.out
~/bin/unphased peblo06.ped -data peblo06.data > peblo06.out
~/bin/unphased peblo07.ped -data peblo07.data > peblo07.out
~/bin/unphased peblo08.ped -data peblo08.data > peblo08.out
~/bin/unphased peblo09.ped -data peblo09.data > peblo09.out
################################################################################
################################################################################
# digit pheno
awk -F","  '{if ($3!="")print $1}' data/BLOCK-DIGIT_23FEB12_4Kuang.csv \
 | grep -v SANGERID > t.ls
# 3131 individuals
grep -f t.ls pe18_12.fam > t.fam
# 1440 individuals

plink --noweb --bfile pe18 --keep t.fam --make-bed --out pedig
t1.py
----------------------------------------
#! /usr/bin/env python3

ifile=open('data/BLOCK-DIGIT_23FEB12_4Kuang.csv')
pheno={}
ifile.readline()
for line in ifile:
    pheno[line.split(',')[0]]=line.split(',')[2]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()
----------------------------------------
t1.py > t.phe

paste pedig.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pedig.fam

plink --noweb --bfile pedig --recode --out pedig

t.head: "T      dig"
awk '{print "M\t" $2}' pedig.map > t.map
cat t.head t.map > pedig.data
~/bin/unphased pedig.ped -data peddig.data > pedel.out
################################################################################

# to get 10 subsets
awk '{print $2}' pedig.bim > snp.ls
# 695464 snps
split -l 69547 snp.ls

plink --noweb --bfile pedig --extract xaa --recode --out pedig00
plink --noweb --bfile pedig --extract xab --recode --out pedig01
plink --noweb --bfile pedig --extract xac --recode --out pedig02
plink --noweb --bfile pedig --extract xad --recode --out pedig03
plink --noweb --bfile pedig --extract xae --recode --out pedig04
plink --noweb --bfile pedig --extract xaf --recode --out pedig05
plink --noweb --bfile pedig --extract xag --recode --out pedig06
plink --noweb --bfile pedig --extract xah --recode --out pedig07
plink --noweb --bfile pedig --extract xai --recode --out pedig08
plink --noweb --bfile pedig --extract xaj --recode --out pedig09

awk '{print "M\t" $2}' pedig00.map > t.map; cat t.head t.map > pedig00.data
awk '{print "M\t" $2}' pedig01.map > t.map; cat t.head t.map > pedig01.data
awk '{print "M\t" $2}' pedig02.map > t.map; cat t.head t.map > pedig02.data
awk '{print "M\t" $2}' pedig03.map > t.map; cat t.head t.map > pedig03.data
awk '{print "M\t" $2}' pedig04.map > t.map; cat t.head t.map > pedig04.data
awk '{print "M\t" $2}' pedig05.map > t.map; cat t.head t.map > pedig05.data
awk '{print "M\t" $2}' pedig06.map > t.map; cat t.head t.map > pedig06.data
awk '{print "M\t" $2}' pedig07.map > t.map; cat t.head t.map > pedig07.data
awk '{print "M\t" $2}' pedig08.map > t.map; cat t.head t.map > pedig08.data
awk '{print "M\t" $2}' pedig09.map > t.map; cat t.head t.map > pedig09.data

~/bin/unphased pedig00.ped -data pedig00.data > pedig00.out
~/bin/unphased pedig01.ped -data pedig01.data > pedig01.out
~/bin/unphased pedig02.ped -data pedig02.data > pedig02.out
~/bin/unphased pedig03.ped -data pedig03.data > pedig03.out
~/bin/unphased pedig04.ped -data pedig04.data > pedig04.out
~/bin/unphased pedig05.ped -data pedig05.data > pedig05.out
~/bin/unphased pedig06.ped -data pedig06.data > pedig06.out
~/bin/unphased pedig07.ped -data pedig07.data > pedig07.out
~/bin/unphased pedig08.ped -data pedig08.data > pedig08.out
~/bin/unphased pedig09.ped -data pedig09.data > pedig09.out
################################################################################




plink --noweb --bfile bed/pewbv --recode --out pewbv
 
awk '{print $2}' bed/pewbv.fam > t.ls
grep -f t.ls pe17_s.evec | awk '{print $2,$3,$4}' > t.cov

paste pewbv.ped t.cov > t.ped; mv t.ped pewbv.ped

awk '{print "M\t" $2}' pewbv.map > t.map
cat t.head t.map t.end > pewbv.data

~/bin/unphased pewbv.ped -data pewbv.data -confounder cov1 cov2 cov3 > pewbv.out



