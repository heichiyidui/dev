################################################################################
# repeat that for the ranking-normaled pheno types                             #
################################################################################
# WBV
awk -F "," '{if ($19!="NA") print $1}' data/Quant_traits_resid.csv | \
 grep -v SANGERID > t.ls
# 1353 sangerids

grep -f t.ls pe18_12.fam > t.fam
# 736 individuals

plink --noweb --bfile pe18 --keep t.fam --make-bed --out pewbv
grep -w -f t.ls pe18_s.evec | awk '{print $2,$3,$4}' > t.cov
#t1.py:
#---------------------------------------
#! /usr/bin/env python3

ifile=open('data/Quant_traits_resid.csv')
pheno={}
ifile.readline()
for line in ifile:
    if line.split(',')[18] != 'NA':
        pheno[line.split(',')[0]]=line.split(',')[18]
ifile.close()

ifile=open('t.fam')
for line in ifile:
    id=line.split()[1]
    print(pheno[id])
ifile.close()
#---------------------------------------
t1.py > t.phe
paste pewbv.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pewbv.fam

plink --noweb --bfile pewbv --recode --out pewbv

 t.head: "T      wbv"
 t.end: 
C       cov1
C       cov2
C       cov3

paste pewbv.ped t.cov > t.ped; mv t.ped pewbv.ped

awk '{print "M\t" $2}' pewbv.map > t.map
cat t.head t.map t.end > pewbv.data

~/bin/unphased pewbv.ped -data pewbv.data -confounder cov1 cov2 cov3 > pewbv.out
################################################################################

# LVV
awk -F "," '{if ($23!="NA") print $1}' data/Quant_traits_resid.csv | \
 grep -v SANGERID > t.ls
# 1345 sangerids
grep -f t.ls pe18_12.fam > t.fam
# 734 individuals

plink --noweb --bfile pe18 --keep t.fam --make-bed --out pelvv
grep -w -f t.ls pe18_s.evec | awk '{print $2,$3,$4}' > t.cov

# change all "[18]" in t1.py to "[22]"
t1.py > t.phe

paste pelvv.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pelvv.fam

plink --noweb --bfile pelvv --recode --out pelvv

t.head: "T      lvv"


paste pelvv.ped t.cov > t.ped; mv t.ped pelvv.ped

awk '{print "M\t" $2}' pelvv.map > t.map
cat t.head t.map t.end > pelvv.data

~/bin/unphased pelvv.ped -data pelvv.data -confounder cov1 cov2 cov3 > pelvv.out

################################################################################
# p3a

awk -F "," '{if ($27!="NA") print $1}' data/Quant_traits_resid.csv | \
 grep -v SANGERID > t.ls
# 1041 sanger ids
grep -f t.ls pe18_12.fam > t.fam
# 501 individuals
plink --noweb --bfile pe18 --keep t.fam --make-bed --out pep3a

grep -w -f t.ls pe18_s.evec | awk '{print $2,$3,$4}' > t.cov

# change all "[16]" in t1.py to "[26]"
t1.py > t.phe

paste pep3a.fam t.phe | awk '{print $1,$2,$3,$4,$5,$7}' > t.out
mv t.out pep3a.fam

plink --noweb --bfile pep3a --recode --out pep3a

t.head: "T      p3a"

paste pep3a.ped t.cov > t.ped; mv t.ped pep3a.ped

awk '{print "M\t" $2}' pep3a.map > t.map
cat t.head t.map t.end > pep3a.data

~/bin/unphased pep3a.ped -data pep3a.data -confounder cov1 cov2 cov3 > pep3a.out

################################################################################
# p3l

awk -F "," '{if ($31!="NA") print $1,$31}' data/Quant_traits_resid.csv | \
grep -w -v "NA" | grep -v SANGERID | awk '{print $1}'> t.ls
# 1053 sanger ids
grep -f t.ls pe18_12.fam > t.fam
# 506 individuals
plink --noweb --bfile pe18 --keep t.fam --make-bed --out pep3l
grep -w -f t.ls pe18_s.evec | awk '{print $2,$3,$4}' > t.cov
# change all "[16]" in t1.py to "[30]"
t1.py > t.phe






