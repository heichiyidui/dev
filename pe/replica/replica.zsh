################################################################################
#             for the replication check vs pgc and sgene                       #
################################################################################

################################################################################
# 1. removing overlapping individuals
#  PGC is removing nothing. 
#  sgene can remove the overlapping vs pgc and us
#  We need to remove the overlapping vs pgc only

#  21 other individuals needed to be removed because of the wrong diagnosis. 
#   (email from Elvira 09 Mar 2012)

awk '{print $1,$2}' ../assoc/pe17.fam | grep -f e21.ls > t.fam
plink --noweb --bfile ../assoc/pe17 --remove t.fam \
 --exclude ../pheno/to_remove.ls --make-bed --out pe19
# 4835 individuals, 695207 SNPs

# t1.py
#--------------------------------------
%#! /usr/bin/env python3
%ifile=open('data/PGC_SCZ_bramon.mesp.overlap.bramon')
%for line in ifile:
%    if line.split()[2].find('_bramon_')==-1:
%        if line.split()[0].find('_bramon_') != -1:
%            print(line.split()[1])
%ifile.close()
#--------------------------------------
# the first column has not our data
t1.py | sort -u > t.out
# 813 sangerids
awk '{print $1,$2}' pe19.fam | grep -f t.out > t.fam
# 717 individuals
plink --noweb --bfile pe19 --remove t.fam --make-bed --out pe19_nopgc
# 4118 individuals, 695207 SNPs


################################################################################
# 2. prepare the ped and data files

awk '{print $2}' pe19.fam > t.ls
grep -f t.ls ../pca/pe17_s.evec | awk '{print $2,$3,$4}' > pe19.cov
awk '{print $2}' pe19_nopgc.fam > t.ls
grep -f t.ls ../pca/pe17_s.evec | awk '{print $2,$3,$4}' > pe19_nopgc.cov

# creat two files: t.head and t.end
# t.head: "A      case"
# t.end:
# C       cov1
# C       cov2
# C       cov3

awk '{print $2}' pe19.bim > t.ls
split -l 6953 t.ls
mv  xaa t00.ls
mv  xab t01.ls
mv  xac t02.ls
...
mv  xdu t98.ls
mv  xdv t99.ls

# t00.sh:
#--------------------------------------
#!/bin/sh
#$ -S /bin/sh
#$ -o  /scratch/home/klinbrc/dev/pe/replica
#$ -e  /scratch/home/klinbrc/dev/pe/replica

cd /home/klinbrc/scratch/dev/pe/replica
 ~/bin/plink --noweb --bfile /home/klinbrc/scratch/dev/pe/replica/pe19 \
 --extract t00.ls  --recode --out t00 

 paste t00.ped pe19.cov > t00.out ; mv t00.out pe19_00.3cov.ped

awk '{print "M\t" $2}' t00.map > t00.out; 
cat t.head t00.out t.end > pe19_00.3cov.data
rm t00.out t00.ped t00.log t00.map
#--------------------------------------
sed 's/t00/t01/g' t00.sh | sed 's/_00/_01/g' > t01.sh
#...
sed 's/t00/t99/g' t00.sh | sed 's/_00/_99/g' > t99.sh

# qsub all t??.sh


# t00.sh:
#--------------------------------------
#!/bin/sh
#$ -S /bin/sh
#$ -o  /scratch/home/klinbrc/dev/pe/replica
#$ -e  /scratch/home/klinbrc/dev/pe/replica

cd /home/klinbrc/scratch/dev/pe/replica
 ~/bin/plink --noweb --bfile /home/klinbrc/scratch/dev/pe/replica/pe19_nopgc \
 --extract t00.ls  --recode --out t00 

 paste t00.ped pe19_nopgc.cov > t00.out ; mv t00.out pe19_nopgc_00.3cov.ped

awk '{print "M\t" $2}' t00.map > t00.out; 
cat t.head t00.out t.end > pe19_nopgc_00.3cov.data
rm t00.out t00.ped t00.log t00.map
#--------------------------------------
sed 's/t00/t01/g' t00.sh | sed 's/_00/_01/g' > t01.sh
#...
sed 's/t00/t99/g' t00.sh | sed 's/_00/_99/g' > t99.sh

################################################################################
# 3. run the analysis

# t00.sh:
#--------------------------------------
#!/bin/sh
#$ -S /bin/sh
#$ -o  /scratch/home/klinbrc/dev/pe/replica
#$ -e  /scratch/home/klinbrc/dev/pe/replica
~/bin/unphased /scratch/home/klinbrc/dev/pe/replica/pe19_00.3cov.ped \
         -data /scratch/home/klinbrc/dev/pe/replica/pe19_00.3cov.data \
         -confounder cov1 cov2 cov3  > pe19_00.out

#--------------------------------------
sed 's/_00/_01/g' t00.sh  > t01.sh
#...
sed 's/_00/_99/g' t00.sh  > t99.sh

# qsub t??.sh

# t00.sh:
#--------------------------------------
#!/bin/sh
#$ -S /bin/sh
#$ -o  /scratch/home/klinbrc/dev/pe/replica
#$ -e  /scratch/home/klinbrc/dev/pe/replica
~/bin/unphased /scratch/home/klinbrc/dev/pe/replica/pe19_nopgc_00.3cov.ped \
         -data /scratch/home/klinbrc/dev/pe/replica/pe19_nopgc_00.3cov.data \
         -confounder cov1 cov2 cov3  > pe19_nopgc_00.out

#--------------------------------------
sed 's/_00/_01/g' t00.sh  > t01.sh
#...
sed 's/_00/_99/g' t00.sh  > t99.sh
# qsub -q long.q,short.q t??.sh

################################################################################
# 4. get the results

# remove the 14 bad snps from the output files

# t.head : CHR SNP BP P
awk '{print $1,$2,$4}' pe19.bim > t.in

grep Likeli res/up_pe19.out | awk '{print $11}' > p.ls
paste t.in p.ls > t.out
cat t.head t.out > up.table

grep Likeli res/up_pe19_nopgc.out | awk '{print $11}' > p.ls
paste t.in p.ls > t.out
cat t.head t.out > up_nopgc.table

#t1.R:
#----------------------------------------
#source("http://www.StephenTurner.us/qqman.r")
#
#results <- read.table("up.table",T);
#png(filename='up_3_man.png',width=1600,height=600);manhattan(results);dev.off();
#png('up_3_qq.png',width=1200,height=1200);gg.qq(results$P);dev.off();
#
#results <- read.table("up_nopgc.table",T);
#png(filename='up_nopgc_3_man.png',width=1600,height=600);manhattan(results);
#dev.off();
#png('up_nopgc_3_qq.png',width=1200,height=1200);gg.qq(results$P);dev.off();
#----------------------------------------

################################################################################
# 5. proxy assoc
awk '{if ($4<0.0001) print $2}' res/up.table > t.ls
# evoker the 57 hits
# all callings are fine

awk '{if ($4<=0.00001) print $2}' res/up.table > t2.ls
# 9 SNPs

# any, grep the proxy SNPs
# t1.py 
#---------------------------------------
#! /usr/bin/env python3

nsnp=57

chr=[]
pos=[]
ifile=open('t.bim')
for line in ifile:
    chr.append(int(line.split()[0]))
    pos.append(int(line.split()[3]))
ifile.close()

ifile=open('pe19.bim')
for line in ifile:
    lchr=int(line.split()[0])
    lpos=int(line.split()[3])
    isProxy=0
    for i in range(nsnp):
        if lchr==chr[i] and lpos < pos[i]+200000 and lpos > pos[i]-200000 : 
            isProxy=1
            continue;
    if isProxy:
        print(line.split()[1])
#---------------------------------------
t1.py > s.ls
plink --noweb --bfile pe19 --extract s.ls --make-bed --out pe19_s

# extract the 4539 of +-200,000 SNPs

awk '{print "plink --noweb --bfile pe19_s --proxy-assoc " $1 " --out " $1}' \
 t.ls > t.sh
source t.sh

cat *.proxy.report > t57_proxy.report
rm *.proxy.report
# pat.ls:
SNP
rs
Proxy

grep -f pat.ls t57_proxy.report > t_proxy.report

awk '{if ($1=="SNP" ||($4<100 && $4>-100 )) print $0}' t_proxy.report | \
awk '{if ($1=="SNP" || $5=="*" || $5>=0.5) print $0}' | \
awk '{if ($8=="P" || $8 < 0.01) print $0}' > t.out

