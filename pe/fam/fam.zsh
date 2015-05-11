################################################################################
# (3) to find the family structure in the PE data                              #
################################################################################

# 1. ibd distances 

# t.ls has the 69802 SNPs from the last SNP pairwise independent selection. 

plink --noweb --bfile pe1 --extract t.ls --maf 0.1 --hwe 1e-3 --make-bed \
 --out pe1_s

# 117 SNPs failed frequency test ( MAF < 0.1 )
# 310 markers to be excluded based on HWE test ( p <= 0.001 )
# 5577 individuals

########################################
# in the brc cluster
split -l 1120 pe1_s.fam

cat xaa xab >  t0_1.fam
cat xaa xac >  t0_2.fam
...
cat xae xad >  t4_3.fam

plink --noweb --bfile pe1_s --keep t0_1.fam --make-bed --out p0_1&
plink --noweb --bfile pe1_s --keep t0_2.fam --make-bed --out p0_2&
...
plink --noweb --bfile pe1_s --keep t4_3.fam --make-bed --out p4_3&

...
~/bin/plink --noweb  --bfile /home/klinbrc/scratch/dev/pe_genome/p4_3 --genome \
 --min 0.05 --out /home/klinbrc/scratch/dev/pe_genome/p4_3

########################################
# given the 20 genome files

cat *.genome | grep -v FID1 | awk '{print $2,$4,$7,$8,$9,$10}' | \
sort -u  > all.genome

# we now have 108278 lines of IID1 IID2 Z0 Z1 Z2 PI_HAT
# 69465 pairs

# The differences between the multiple entries of the same pair are tiny.
# to remove duplicates

t1.py :
----------------------------------------
#! /usr/bin/env python3

ifile=open('all.genome')

comIds=[]
line = ifile.readline()
if line.split()[0]>line.split()[1] :
    comIds.append(line.split()[0]+line.split()[1])
else:
    comIds.append(line.split()[1]+line.split()[0])
print(line[:-1])

for line in ifile:
    if line.split()[0]>line.split()[1] :
        comId=line.split()[0]+line.split()[1]
    else:
        comId=line.split()[1]+line.split()[0]
    if comId in comIds:
        continue;
    
    print(line[:-1])
    comIds.append(comId)

----------------------------------------

t1.py | sort -g -r +5 > t.out; cp t.out all.genome
# ibd distancec of 69465 pairs of individuals

################################################################################
# 2. remove the duplication of individuals 
# (some might be identical twins)

# first remove the MT,X,Y and XY chrs
awk '{if ($1>22) print $2}' ../../PE/calls/pe1.bim > t.ls

plink --noweb --bfile ../../PE/calls/pe1 --exclude t.ls --make-bed --out pe2
# removed 38895 SNPs
# left 890661 SNPs, 5577 individuals

# checking the missingness of individuals

plink --noweb --bfile pe2 --missing --out pe2_missing
tail +2 pe2_missing.imiss | awk '{print $2,$6}' > t.missing

# 87 pairs of duplicates
sort -g -r +4 all.genome | head -n 87 | awk '{print $1,$2}' > t.in
# but only 157 individules
awk '{print $1,"\n",$2}' t.in | sed 's/ //g' | sort -u | wc 

# list the individuals of the duplicates pairs with higher missingness
t1.py
----------------------------------------
#! /usr/bin/env python3

mfile=open('t.missing')
missingness={}
for line in mfile:
    missingness[line.split()[0]]=float(line.split()[1])
mfile.close()

ifile=open('t.in')
for line in ifile:
    if missingness[line.split()[0]]>missingness[line.split()[1]]:
        print(line.split()[0])
    else:
        print(line.split()[1])

----------------------------------------

t1.py | sort -u > dup_to_remove.ls

# OK, 81 individuals to be removed.

grep -w -f dup_to_remove.ls pe2.fam > t.fam 

plink --noweb --bfile pe2 --remove t.fam --make-bed --out pe3
# 5496 individuals,  890661 SNPs

grep -w -v -f dup_to_remove.ls all.genome > t.out; mv t.out all.genome
# 68338 pairs of relatives

# about the 40 heidenbuge
# get the 40 sangerids into hei40.ls
grep -w -f hei40.ls ../../PE/calls/pe1.fam | awk '{print $2}' > hei36.ls
# 36 in the pe1 set
grep -w -f hei36.ls t.in
#WTCCCT628592 WTCCCT628591
#WTCCCT628639 WTCCCT628638
#WTCCCT628582 WTCCCT628583
#WTCCCT628585 WTCCCT628584
#WTCCCT628573 WTCCCT628574
#WTCCCT628554 WTCCCT628555
#WTCCCT628649 WTCCCT628648
#WTCCCT628562 WTCCCT628563
#WTCCCT628628 WTCCCT628627
#WTCCCT628593 WTCCCT628594
#WTCCCT628580 WTCCCT628581
# 11 pairs
grep -w -f hei36.ls pe3.fam| awk '{print $2}' > hei25.ls
# 25 in the pe3 set

################################################################################
# 3. correct the phenotypes, 

# use the Phenotype 4 Kuang 8797 - 18Nov11.xls
# put the SangerId and the Phenotype columns into t.in
# the 40 'Exclude' are now all controls

t1.py
----------------------------------------
#! /usr/bin/env python3

index=open('t.in')
caseInd=[]
contInd=[]
for line in index:
    if line.split()[1]=='Patient':
        caseInd.append(line.split()[0])
    if line.split()[1]=='Control':
        contInd.append(line.split()[0])

ifile=open('pe3.fam')
for line in ifile:
    col=line.split()
    if col[1] in caseInd:
        print(col[0],col[1],col[2],col[3],col[4],'2')
    elif col[1] in contInd:
        print(col[0],col[1],col[2],col[3],col[4],'1')
    else: # relatives, exclude etc
        print(col[0],col[1],col[2],col[3],col[4],'-9')
----------------------------------------
t1.py > t.out ; mv t.out pe3.fam

# 1035 relatives, 2967 controls and 1494 cases.

################################################################################
# 4. check relatives
# we have the gender information for all individuals

# we know for different relationships:
#                 z1     z2     z3
# identical        0      0      1
# parent-child     0      1      0
# siblings         0.25   0.5    0.25 
# half-siblings    0.5    0.5    0
# uncle-niece      0.5    0.5    0
# First cousins    0.75   0.25   0

# cut z2 at 0.75, the first 1242 are considered as parent/offspring pairs
sort -g -r +3 all.genome | head -n 1242 | awk '{print $1,$2}' > po.pairs

awk '{print $1,"\n",$2}' po.pairs | sed 's/ //g' | sort -u | wc
#1641 individuals

# cut parent-offspring pairs, 
# set z2 and z3 thresholds at 0.375 and 0.125

sort -g -r +3 all.genome | tail +1243 | \
awk '{if ($4>0.375 && $5>0.125) print $1,$2}' > si.pairs
# 756 pairs

awk '{print $1,"\n",$2}' si.pairs | sed 's/ //g' | sort -u | wc
# 1160 individuals

# use the MRI data 4 Kuang at 21Nov 11_PAPER 8797.xls table
# get the id and age columns into t.in
# update using the DOB GROUP from Ruud final 28Nov11.xls table
# got the gender information from pe3.fam

# 553 individuals have no age info

awk '{print $1,"\n",$2}' po.pairs | sed 's/ //g' | sort -u > t.ls
grep -w -f t.ls t.in | sort -u | awk '{if ($3==-9) print $1}' > to_be_aged.ls

# 15 individuals of parent/offspring relationships need ages to be assigned.

################################################################################
# 5. assign ages

# to check and correct available ages

# found 2.0 3.0 4.0 6.0 8.0 years old fathers
# and 11.0 12.0 13.0 13.0 years old mothers, and maybe the 48 years old one...
# these should be fixed

# on average, mothers should be 28 years older, fathers 31

# 5.1 the parents should be at least 13 years elder
t1.py
----------------------------------------
#! /usr/bin/env python3
import math

index=open('t.in')
sex={}
age={}
for line in index:
    col=line.split()
    sex[col[0]]=int(col[1])
    age[col[0]]=float(col[2])
index.close()

ifile=open('po.pairs')
for line in ifile:
    col=line.split()
    if age[col[1]] == -9 or age[col[0]] == -9:
        continue
    agediff=abs( age[col[1]] - age[col[0]])
    if agediff < 14:
        print(col[0],sex[col[0]],age[col[0]])
        print(col[1],sex[col[1]],age[col[1]],'\n')
ifile.close()
----------------------------------------
# family 1
 parents
 WTCCCT535377 1 57.00
 WTCCCT534930 2 53.0 
 siblings
 WTCCCT535120 1 41.0
 WTCCCT535327 2 30.00
 fix (using the age of the sibling)
 WTCCCT535120 -> 30
# family 2
 grandfather
 WTCCCT526544 1 61.0 
 mother
 WTCCCT526543 2 53.0
 daughter
 WTCCCT526542 2 21.00
 fix (using the average age of the father and the daughter)
 WTCCCT526543 -> 41
# family 3
 parents
 WTCCCT525316 1 22
 WTCCCT525317 2 46.0
 siblings
 WTCCCT525315 1 20.0
 WTCCCT525318 1 49.0 
 fix (switch the ages of the father and a son)
 WTCCCT525316 -> 49
 WTCCCT525318 -> 22
# family 4
 mother
 WTCCCT535386 2 63.0 
 siblings
 CCC2_PE815151 2 50.0
 WTCCCT535217  1 32.00
 CCC2_PE814886 1 36.00
 fix (using the average age of the siblings)
 CCC2_PE815151 -> 34
# family 5
 parents
 WTCCCT526949 1 36.0
 WTCCCT526937 2 63.00
 sisters
 WTCCCT525134 2 32.00
 WTCCCT525091 2 30.0 
 fix (using the age of the mother)
 WTCCCT526949 -> 63
# In family 6, 7 mother/offspring age differences are 11 and 13. No other info. 
# Didn't change.
# 6 assigned

########################################
# 5.2 assign the missing ages

awk '{print $1,"\n",$2}' po.pairs | sed 's/ //g' | sort -u > t.ls
grep -w -f t.ls t.in | sort -u | awk '{if ($3==-9) print $1}' > to_be_aged.ls
# 15 individuals with parent/offspring relationships 
# need to be assigned with ages.

# to find the parent/offspring and siblings:
t1.py
----------------------------------------
#! /usr/bin/env python3

index=open('t.in')
sex={}
age={}
for line in index:
    col=line.split()
    sex[col[0]]=int(col[1])
    age[col[0]]=float(col[2])
index.close()

ifile=open('to_be_aged.ls')
for line in ifile:
    id=line.strip()
    print('\nmissing:',id)
    
    pfile=open('po.pairs')
    for line in pfile:
        col=line.split()
        pid=''
        if col[0]==id:
            pid=col[1]
        if col[1]==id:
            pid=col[0]
        if pid:
            print('po',pid,sex[pid],age[pid])
    pfile.close()
    
    sfile=open('si.pairs')
    for line in sfile:
        col=line.split()
        sid=''
        if col[0]==id:
            sid=col[1]
        if col[1]==id:
            sid=col[0]
        if sid:
            print('si',sid,sex[sid],age[sid])
    sfile.close()
ifile.close()
----------------------------------------

# The average age is 43.
# The average age difference between parent/child pair is 29.

# The rules are simple. 
# The age of some siblings are known, use the average age of known siblings.

missing: WTCCCT534986
po WTCCCT534914 1 33.0
si CCC2_PE815166 2 58.0
si WTCCCT534951 2 49.0

fix: WTCCCT534986 -> 54

missing: WTCCCT535291
po WTCCCT535626 2 61.0
si WTCCCT535197 2 38.0
si WTCCCT535578 1 36.0
si WTCCCT535613 1 39.0

fix: WTCCCT535291 -> 37

missing: WTCCCT535339
po WTCCCT536231 1 61.0
si WTCCCT535532 1 39.0

fix: WTCCCT535339 -> 39

3 assignments

################################################################################
# 6. get the new pedegree files

# 6.1 assign parents
t1.py 
----------------------------------------
#! /usr/bin/env python3

index=open('t.in')
sex={}; age={}
father={}; mother={}
for line in index:
    col=line.split()
    sex[col[0]]=int(col[1])
    age[col[0]]=float(col[2])
    father[col[0]]='0'
    mother[col[0]]='0'
index.close()

ifile=open('po.pairs')
for line in ifile:
    col=line.split()
    if age[col[1]] == -9 or age[col[0]] == -9:
        continue
    idp=col[1] # id of parent
    ido=col[0] # id of offspring
    if age[idp] < age[ido]:
        idp=col[0]
        ido=col[1]
    if sex[idp]==1:
        father[ido]=idp
    else:
        mother[ido]=idp
ifile.close()

ifile=open('pe3.fam')
for line in ifile:
    col=line.split()
    col[2]=father[col[1]]
    col[3]=mother[col[1]]
    print(col[0],col[1],col[2],col[3],col[4],col[5])
ifile.close()
----------------------------------------
t1.py > t.out; cp t.out pe3.fam

########################################
# 6.2 correct family ids

awk '{print "F" NR, $2,$3,$4,$5,$6}' pe3.fam > t.out
cp t.out pe3.fam

t1.py 
----------------------------------------
#! /usr/bin/env python3

ifile=open('pe3.fam')
fam={} # family ids, without the leading 'F', as numbers
for line in ifile:
    fam[line.split()[1]]=int(line.split()[0][1:])
ifile.close()

for i in range(10):
    ifile=open('pe3.fam')
    for line in ifile:
        col=line.split()
        id=col[1]
        father=col[2]; mother=col[3]

        if father != '0':
            if fam[father] < fam[id]:
                fam[id]=fam[father]
            if fam[id] < fam[father]:
                fam[father]=fam[id]
        if mother !='0':
            if fam[mother] < fam[id]:
                fam[id]=fam[mother]
            if fam[id] < fam[mother]:
                fam[mother]=fam[id]
    ifile.close()

ifile=open('pe3.fam')
for line in ifile:
    col=line.split()
    id=col[1];
    print('F'+str(fam[id]),id,col[2],col[3],col[4],col[5])
ifile.close()
----------------------------------------
t1.py > t.out; cp t.out pe3.fam

awk '{print $1}' pe3.fam | sort | uniq -d | wc
# got 531 families

awk '{print $1}' pe3.fam |sort| uniq -c | awk '{print $1}' | sort -g -r > t.out
# of all 4402 families (including singleton)
# 5 have 6 members, 29 have 5 members, 117 with 4, 222 with 3, 158 with 2
# 3871 singletons

########################################
# 6.3 for pedstats, we need both parents to be present, or none at all
cp pe3.fam pe3_full.fam
awk '{$1="F" NR;if ($3==0)$4=0;if($4==0)$3=0; print $0}' pe3.fam > pe3_part.fam

# re-assign family numbers
sed -i 's/pe3.fam/pe3_part.fam/g' t1.py
t1.py > t.out; cp t.out pe3_part.fam

awk '{print $1}' pe3_part.fam | sort | uniq -d | wc

awk '{print $1}' pe3_part.fam |sort| uniq -c | \
 awk '{print $1}' | sort -g -r > t.out

grep 6 t.out | wc 
grep 5 t.out | wc
grep 4 t.out | wc
grep 3 t.out | wc
grep 1 t.out | wc

# 225 families
# 5 with 6 members, 19 with 5 members, 89 with 4 members, 115 with 3 members
# 4682 singletons
