################################################################################
#                    To impute the QCed PE GWAS data                           #
################################################################################

################################################################################
# 1. To prepare pe19 GEN files

########################################
# 1.1 to update map to build 37 

mkdir bed 
cp ../pe_pheno/data/pe19.??? bed

# 4835 individuals
# 694673 SNPs
# AT+TA+CG+GC = 21708+22265+32420+31327 = 107720
# 15.5 percent

wget www.well.ox.ac.uk/~wrayner/strand/GenomeWideSNP_6.na32-b37.strand.zip
unzip GenomeWideSNP_6.na32-b37.strand.zip

awk '{if ($5=="---") print $2}' GenomeWideSNP_6.na32.strand.txt > t1.ls 
# 3339 SNPs have no direction, of which 1846 SNPs have no chr and position info.

tail -n +2 GenomeWideSNP_6.na32.strand.txt | awk '{print $2}' > t.ls
grab -f t.ls -v -c 2 bed/pe19.bim | awk '{print $2}' > t2.ls
# 464 SNPs in pe19 but not in the strand file 

awk '{if ($3=="X" || $3=="Y" || $3=="MT") print $2}' \
    GenomeWideSNP_6.na32.strand.txt > t3.ls 
# 39173 SNPs in X or Y or MT 

# remove them 
cat t1.ls t2.ls t3.ls > t.ls
plink --noweb --bfile bed/pe19 --exclude t.ls --make-bed --out pe19 
# 694673 to 693530 SNPs, 1143 SNPs removed

# update SNP chr and position 
t1.py 
#--------------------------------------
#!/usr/bin/env python3 

snp_chr={}
snp_pos={}
ifile=open('GenomeWideSNP_6.na32.strand.txt')
ifile.readline()
for line in ifile:
    cols=line.split()
    snp_chr[cols[1]]=cols[2]
    snp_pos[cols[1]]=cols[3]
ifile.close()

ifile=open('pe19.bim')
for line in ifile:
    cols=line.split()
    snp_id=cols[1]
    cols[0]=snp_chr[snp_id]
    cols[3]=snp_pos[snp_id]
    print('\t'.join(cols))
ifile.close()
#--------------------------------------

t1.py > t.out
mv t.out pe19.bim 

# because of the implicit order changed from re-mapping 
# need to rewrite the files
plink --noweb --bfile pe19 --make-bed --out pe19

# check SNPs sorted by position
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    chr_int=$(echo $chr | sed 's/^[0]//')
    awk -v chr_int=$chr_int '{if ($1==chr_int) print $4}' pe19.bim > t.ls 
    sort -g t.ls > t2.ls 
    echo $chr 
    wc -l t.ls 
    diff t.ls t2.ls | wc
    uniq -d t2.ls | wc 
done 

#######################################
# 1.2 flip the strand of snps
tail -n +2 GenomeWideSNP_6.na32.strand.txt | awk '{if ($5=="-") print $2}' >t.ls

plink --noweb --bfile pe19 --flip t.ls --make-bed --out pe19 
# Flipped strand of 348688 SNPs

# checked the allele frequencies of 10 SNPs vs dbSNP. All fine. 

########################################
# 1.3 convert to gen files

# Now not neccessary with shapeit 

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N gen_convert
#$ -cwd
#$ -q long.q,bignode.q,short.q

chr=$1                                   # 01 02 ... 22
chr_int=$(echo $chr | sed 's/^[0]//')    #  1  2 ... 22

~/bin/plink --noweb --bfile pe19 --chr $chr_int --recode   --out pe19_$chr 
~/bin/plink --noweb --bfile pe19 --chr $chr_int --make-bed --out pe19_$chr 

/share/bin/gtool -P --ped pe19_${chr}.ped --map pe19_${chr}.map \
  --og pe19_${chr}.gen --os pe19_${chr}.sample 

#--------------------------------------
# end of t.sh 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    qsub t.sh  $chr 
done 

########################################
# 1.4 re-check strand with shapeit
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    
    chr_int=$(echo $chr | sed 's/^[0]//')

ref_dir=/home/klinbrc/scratch/ref/1000GP_Phase3
ref_head=1000GP_Phase3
~/bin/shapeit.v2.r790.linux.x64 \
        -check \
        -B pe19_${chr} \
        -M  $ref_dir/genetic_map_chr${chr_int}_combined_b37.txt \
         --input-ref $ref_dir/${ref_head}_chr${chr_int}.hap.gz \
                $ref_dir/${ref_head}_chr${chr_int}.legend.gz \
                $ref_dir/${ref_head}.sample \
        --output-log pe19_${chr}.alignments
done

cat pe19_??.alignments.snp.strand | grep Missing | wc
# 2278 SNPs in pe19, but not in 1000G panel. Remove them.
cat pe19_??.alignments.snp.strand | grep Strand | wc
# 1521 SNPs to flip

rm t.ls
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    awk '{if ($1=="Missing") print $4}' pe19_$chr.alignments.snp.strand >> t.ls
done

plink --noweb --bfile pe19 --exclude t.ls --make-bed --out pe19

rm t.ls
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    awk '{if ($1=="Strand") print $4}' pe19_$chr.alignments.snp.strand >> t.ls
done

plink --noweb --bfile pe19 --flip t.ls --make-bed --out pe19
# 691252 SNPs left

# redo the conversion then

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    qsub t.sh  $chr 
done 

# Can still find some miss-aligned SNPs now with shapeit. 
# Don't bother. 

################################################################################
# 2. Phasing for impute

#######################################
# 2.2 Phasing

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N shapeit
#$ -cwd
#$ -q long.q,bignode.q,short.q

#$ -pe multi_thread 6
#$ -l h_vmem=6G

chr=$1                                   # 01 02 ... 22
chr_int=$(echo $chr | sed 's/^[0]//')    #  1  2 ... 22

ref_dir=/home/klinbrc/scratch/ref/1000GP_Phase3
~/bin/shapeit.v2.r790.linux.x64 \
    --thread 6 \
    --input-bed pe19_${chr}.bed pe19_${chr}.bim pe19_${chr}.fam \
    --input-map $ref_dir/genetic_map_chr${chr_int}_combined_b37.txt \
    --output-max pe19_${chr}.phased.haps pe19_${chr}.phased.sample \
    --duohmm
#---------------------------------------
# end of t.sh 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    qsub t.sh $chr 
done 

# 2.5 hours for chr 22 
# 24  hours for chr 1

########################################
# 2.2 to split maps to chunks and pieces

t1.py 
#--------------------------------------
#!/usr/bin/env python3

import sys

PIECES_SNP_SIZE  = 1200     # roughly that many SNPs in any pieces
PIECES_LENGTH    = 4500000  # rough piece size in bp
BUFFER_SIZE      = 250000   # default buffer size
MIN_PIECE_SIZE   = 50       # in a chunk, should have at least 51 SNPs 
#######################################
# to read SNP positions and chr 

ifile=open(sys.argv[1]) # input bim or map file
#ifile=open('pe19_22.bim')
snp_pos=[]  # SNP bp position

line=ifile.readline()
cols=line.split()
chr='{:02d}'.format(int(cols[0])) # 01 02 ... 22

snp_pos.append(int(cols[3]))

for line in ifile:
    snp_pos.append(int(line.split()[3]))
ifile.close()

#######################################
# to cut chromsomes into chunks
# SNPs in a chunk can not be seperated by more than PIECES_LENGTH / 4

chunks=[]
chunk=[snp_pos[0]]
for snp_po in snp_pos[1:]:
    if snp_po - chunk[-1] > PIECES_LENGTH / 4:
        chunks.append(chunk)
        chunk=[snp_po]
    else:
        chunk.append(snp_po)
chunks.append(chunk)

# print('\nchr:',chr, '\tnum_of_chunks:',len(chunks))

#######################################
# to cut chunks into pieces

pieces=[]
for chunk in chunks:
    # print('\tchunk_begin:',chunk[0],\
    #       '\tchunk_end:',  chunk[-1],\
    #       '\tchunk_length:',chunk[-1]-chunk[0]+1,\
    #       '\tSNPs_in_chunk:',len(chunk))
    if len(chunk) < MIN_PIECE_SIZE:
        continue
    len_chunk=chunk[-1]-chunk[0]
    num_pieces=len_chunk//PIECES_LENGTH + 1
    # print('\t\tnum_of_pieces_in_chunk:',num_pieces,\
    #       '\taverage_SNPs_in_pieces:',len(chunk)//num_pieces)
    
    size_pieces=len(chunk)//num_pieces
    for i in range(0,len(chunk),size_pieces):
        pieces.append(chunk[i:i+size_pieces])

for i in range(7):
    new_pieces=[]
    for piece in pieces:
        if len(piece) < MIN_PIECE_SIZE:
            continue
        if piece[-1] - piece[0] < PIECES_LENGTH:
            new_pieces.append(piece)
            continue
        
        dist=[] # distances between snps
        for j in range(len(piece)-1):
            dist.append(piece[j+1]-piece[j])
        sep=dist.index(max(dist))# the seperation point
        
        if sep < len(piece)*0.1 or sep > len(piece)*0.9:
            sep=len(piece)//2
        
        new_pieces.append(piece[0:sep+1])
        new_pieces.append(piece[sep+2:])
    pieces=new_pieces

#######################################
# to get the begins and ends of pieces

begins=[]
ends=[]
begins.append(max(pieces[0][0]-BUFFER_SIZE,1))

for i in range(len(pieces)-1):
    gap_begin = pieces[i][-1]
    gap_end   = pieces[i+1][0]
    
    ends.append  (min( (gap_begin+gap_end)//2  ,gap_begin+BUFFER_SIZE))
    begins.append(max( (gap_begin+gap_end)//2+1,gap_end  -BUFFER_SIZE))
ends.append(pieces[-1][-1]+BUFFER_SIZE)

#######################################
# to output the command

for i in range(len(begins)):
    piece_code='{:03d}'.format(i)
    print('qsub impu.sh',chr,piece_code,begins[i],ends[i])

#--------------------------------------
# end of t1.py 

rm t.out 
for ifile in pe19_??.bim ; do 
    t1.py $ifile >> t.out
done

#
# 866 pieces 
# length ranges from 0.5M to 4.9M 

#######################################
# 2.3 imputation 
# reference from 

# https://mathgen.stats.ox.ac.uk/\
#impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html
# October 2014

impu.sh
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N impute2
#$ -cwd
#$ -q long.q,bignode.q
#$ -l h_vmem=10G

chr=$1
chunk_code=$2
chunk_begin=$3
chunk_end=$4

chr_int=$(echo $chr | sed 's/^[0]//')
ref_dir=/home/klinbrc/scratch/ref/1000GP_Phase3
ref_head=1000GP_Phase3

/share/bin/impute2 \
    -use_prephased_g \
    -known_haps_g pe19_${chr}.phased.haps  \
    -m $ref_dir/genetic_map_chr$chr_int\_combined_b37.txt \
    -h $ref_dir/${ref_head}_chr${chr_int}.hap.gz \
    -l $ref_dir/${ref_head}_chr${chr_int}.legend.gz \
    -int $chunk_begin $chunk_end \
    -Ne 20000 \
    -o pe20_$chr\_$chunk_code.imputed 

#----------------------------------------
# end of impu.sh 

source t.out 

################################################################################
# 3 QC of the imputed data

########################################
# 3.1 check and remove duplications in each chr

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    head -n 1 pe20_${chr}_001.imputed_info > t.head 
    cat pe20_${chr}_*.imputed_info | grep -v snp_id > t.out 
    cat t.head t.out > pe20_${chr}.imputed_info 
    tail -n +2 pe20_${chr}.imputed_info | awk '{print $2}' \
        | sort | uniq -d > dup_${chr}.ls 
    wc dup_${chr}.ls 
done 

# chr 18 19 21 22 have no duplications
# mostly only a few duplications per chr 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 20 ; do  
    grab -v -f dup_${chr}.ls -c 2 pe20_${chr}.imputed_info > t.out 
    mv t.out pe20_${chr}.imputed_info
done
# 14114556 SNPs 

# replace the leading '---' with chr
# note the using of double quotes here. 
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    chr_int=$(echo $chr | sed 's/^[0]//')
    sed -i "s/\-\-\-/$chr_int/" pe20_$chr.imputed_info
done 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    chr_int=$(echo $chr | sed 's/^[0]//')
    cat pe20_$chr\_*.imputed | \
        grep -v -w -f dup_$chr.ls | \
        sed "s/\-\-\-/$chr_int/" > t_$chr.gen
done 

#######################################
# 3.2 basic filter

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    qctool -g t_$chr.gen -maf 0.001 0.5 -info 0.3 1 -og pe20_$chr.bgen 
done 

# do missing-call-rate here
# reuse the qt.sh

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    /share/bin/qctool -g pe20_$chr.bgen -snp-stats pe20_${chr}.stats
done 


#######################################
# 3.3 check allele frequency

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    /share/bin/qctool -g pe20_$chr.bgen -snp-stats pe20_${chr}.stats
done 

ref_dir=/scratch/data/1000genomes/ALL_1000G_phase1integrated_v3_impute_macGT1
legend_head=ALL_1000G_phase1integrated_v3_chr

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
        chr_int=$(echo $chr | sed 's/^[0]//');
        cp $ref_dir/$legend_head$chr_int\_impute_macGT1.legend.gz .
        gunzip $legend_head$chr_int\_impute_macGT1.legend.gz  
        mv  $legend_head$chr_int\_impute_macGT1.legend ref_$chr.legend
done

cat ref_??.legend | grep -v eur.aaf | \
    awk '{if (rand()<0.001) print $1,$3,$4,$8}' > ref.freq

awk '{print $1}' ref.freq > t.ls

rm t.out
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    grab -f t.ls -c 2 pe20_$chr.stats >> t.out
done

awk '{print $2,$6,$5, ($9+$10*0.5)/($9+$10+$11)}' t.out > imp.freq

awk '{print $1}' imp.freq > t.ls
grab -f t.ls ref.freq > t.out
sort_table -f t.ls t.out > ref.freq

paste ref.freq imp.freq > t.in
# 14307 SNPs 

rm t.out 
awk '{if ($2==$6 && $3==$7) print $4,$8}' t.in > t.out
awk '{if ($2==$7 && $3==$6) print $4,1-$8}' t.in > t.out

# 14307 SNPs 

awk '{print $1-$2}' t.out > t.ls 
listdis t.ls > t.dat

# mean 0.0001 std 0.0139

################################################################################
# 4. association study

#######################################
# 4.1 get the sample file 
awk '{print $2,$3,$4,$5}' bed/pe19_3.cov  > t1.in

tail -n +3 pe19_22.sample > t2.in

sort_table -f t.ls t1.in > t3.in
paste t2.in t3.in | awk '{print $1,$2,$3,$7,$8,$9,$4,$5}' > t.out

awk '{if ($8 != -9) $8=$8-1; print $0}' t.out > t.in

t.head 
#--------------------------------------
ID_1 ID_2 Missing cov1 cov2 cov3 gender status
0 0 0 C C C D P
#--------------------------------------
# end of t.head 

cat t.head t.in > pe20.sample 

#######################################
# 4.2 SNPTEST

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -N snptest
#$ -cwd
#$ -q long.q,bignode.q

chr=$1

/share/bin/snptest \
 -data pe20_$chr.bgen pe20.sample \
 -o t_$chr.out \
 -frequentist 1 \
 -method score \
 -pheno status
 
#--------------------------------------
# end of t.sh 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    qsub t.sh $chr 
done

t.head 
#--------------------------------------
SNP CHR BP P
#--------------------------------------
# end of t.head

cat t.head > assoc.p_table 
for ifile in t_??.out ; do 
    tail -n +2 $ifile | awk '{print $2,$1,$4,$20}' >> assoc.p_table 
done 
# 14114556 SNPs

#######################################
# 4.3 SNP QC filtering 

# maf > 0.01
# info >  0.4
# posterior > 0.9

rm t.ls
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do  
    awk '{if ($4<0.99 && $4>0.01 && $5>0.4 && $6>0.9) print $2}' \
        info/pe20_$chr.imputed_info  >> t.ls
done

head -n 1 assoc.p_table > sub_assoc.p_table
grab -f t.ls assoc.p_table >> sub_assoc.p_table
# 8864591 SNPs after filtering

# after info > 0.9 filter, 13799347 SNPs

awk '{if ($4<0.000001) print $1}' assoc.p_table > t.ls 
# 5904 with P values less than 1^-06
awk '{if ($4<0.000001) print $1}' sub_assoc.p_table > t.ls 
# 268 after filtering

awk '{if ($4<0.00001 ) print $1}' assoc.p_table > t.ls 
# 11996 with P values less than 1^-05
awk '{if ($4<0.00001 ) print $1}' sub_assoc.p_table > t.ls 
# 1300 after filtering

#######################################
# 4.4 QQ and Manhattan plots

R
#--------------------------------------
install.packages('qqman')
#--------------------------------------

R
#--------------------------------------
library(qqman)
results=read.table('sub_01_assoc.p_table',header=TRUE)

png('t_01_m.png', width=1400)
manhattan(results)
dev.off()

png('t_01_qq.png')
qq(results$P)
dev.off()

results=read.table('sub_02_assoc.p_table',header=TRUE)

png('t_02_m.png', width=1400)
manhattan(results)
dev.off()

png('t_02_qq.png')
qq(results$P)
dev.off()


results=read.table('sub_03_assoc.p_table',header=TRUE)

png('t_03_m.png', width=1400)
manhattan(results)
dev.off()

png('t_03_qq.png')
qq(results$P)
dev.off()
#--------------------------------------

# chr6:384955:I 6 384955 5.0402e-08 looks good 


# the top hit rs72817292 with maf 0.011 is bad, remove it 
# re-plot
# QQ plot shows it's totally inflated. 

# maf > 0.01 , > 0.02 and > 0.03
# 8864591  35458364 274391957 sub_01_assoc.p_table
# 7802704  31210816 240822816 sub_02_assoc.p_table
# 7097945  28391780 218647927 sub_03_assoc.p_table

# THE QQ plots are always inflated...




#######################################
# gtool conversion to plink files

#$ -l h_vmem=6G

chr=$1
chri=$2

sed -i 's/NA //' pe23_${chr}.gen 

~/bin/gtool -G --g pe23_${chr}.gen --s pe19.sample \
 --ped pe23_${chr}.ped --map pe23_${chr}.map \
 --phenotype Phenotype --sex Gender \
 --threshold 0.9 --chr $chri \
 --snp --log gtlog_$chr 

sed -i 's/\tN N/\t0 0/g' pe23_${chr}.ped

plink --noweb --file pe23_${chr} --make-bed --out pe24_${chr} 
#-------------
qsub qt.sh 01 1
qsub qt.sh 02 2
#...
qsub qt.sh 22 22

for fam in pe24_??.fam; do
    cp bed/pe19.fam $fam
done

################################################################################
# 5. plink assoc

chr=$1

~/bin/plink --noweb --bfile pe24_${chr} --covar bed/pe19_3.cov \
   --linear --standard-beta --ci 0.95 --hide-covar \
   --out pe24_${chr} 

################################################################################
# 6. qq and man plots
awk -v OFS='\t' '{print $1,$2,$3,$4,$7,$8,$12}' pe24_??.assoc.logistic > t.out

head -n 1 t.out > t.head
grep -v CHR t.out | grep -v MERGED_DEL > t.end        
cat t.head t.end > pe24_log.assoc
R --vanilla --slave --args table=pe24_log.assoc fig=pe24 < ~/bin/qqman.r 

################################################################################
# the end                                                                      #
################################################################################
