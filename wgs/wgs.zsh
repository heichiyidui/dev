################################################################################
#               NEXT GENERATION WHOLE GENOME SEQUENCING SCRIPTS                #
################################################################################

################################################################################
# 1. Input                                                                     #
################################################################################

# low coverage sequencing data from Illumina
# 8 samples: MF_1348 MF_2874 MF_3925 MF_3934 MF_4057 MF_4300 MF_4306 MF_5269
# all Schiz cases 
# 2 cells: FlowCell_D0L9C and FlowCell_D0LAC
# 8 x 2 lanes 
# ~50 R1 and R2 fastq.gz file pairs each lane 
# 34 ~ 52 G files each lane 
# 664G altogether , ~330G per cell 
# originaly in /scratch/project/grouppowell/wholegenome
# latter archived into /archive/klinbrc/FlowCell_D0L[9A]C

########################################
# 1.1 to combine multiple fastq files of the same lane 

cd /home/klinbrc/scratch/dev/wgs
mkdir fastq 

/archive/klinbrc/FlowCell_D0L9C/Sample_MF_1348/*R1*.gz
/archive/klinbrc/FlowCell_D0L9C/Sample_MF_1348/*R2*.gz

for cell in D0L9C D0LAC; do 
    for sample in 1348 2874 3925 3934 4057 4300 4306 5269; do 
        gunzip /archive/klinbrc/FlowCell_${cell}/Sample_MF_${sample}/*.gz 
    done
done

for cell in D0L9C D0LAC; do 
    for sample in 1348 2874 3925 3934 4057 4300 4306 5269; do 
        qsub t1.sh $cell $sample 
    done
done

for cell in D0L9C D0LAC; do 
    for sample in 1348 2874 3925 3934 4057 4300 4306 5269; do 
cat /archive/klinbrc/FlowCell_${cell}/Sample_MF_${sample}/*_R1_*.fastq \
    > fastq/${cell}_MF_${sample}_R1.fastq 
cat /archive/klinbrc/FlowCell_${cell}/Sample_MF_${sample}/*_R2_*.fastq \
    > fastq/${cell}_MF_${sample}_R2.fastq 
    done
done

# 16 pairs of R1 R2 fastq files 
# 1.8 T 

################################################################################
# 2. QC and novoalign                                                          #
################################################################################

# 2.1 fastqc

for fqfile in fastq/*.fastq; do fastqc $fqfile; done &
mkdir fastqc
fastq/*fastqc fastqc

grep WARN fastqc/*fastqc/summary.txt 
# all fastq files got Per sequence GC content and Kmer Content warnings

grep -v PASS fastqc/*fastqc/summary.txt 
# no other error messages

for fqfile in FlowCell_D0LAC/*/*.fastq; do fastqc $fqfile; done &

# 39 FAIL Per base sequence quality in D0L9C, all R2
# MF_1348_NoIndex_L008_R2_029.fastq
# MF_1348_NoIndex_L008_R2_036.fastq
# MF_1348_NoIndex_L008_R2_043.fastq
# MF_1348_NoIndex_L008_R2_044.fastq
# MF_2874_NoIndex_L001_R2_016.fastq
# MF_2874_NoIndex_L001_R2_017.fastq
# MF_2874_NoIndex_L001_R2_025.fastq
# MF_2874_NoIndex_L001_R2_033.fastq
# MF_2874_NoIndex_L001_R2_041.fastq
# MF_2874_NoIndex_L001_R2_042.fastq
# MF_2874_NoIndex_L001_R2_050.fastq
# MF_3925_NoIndex_L003_R2_047.fastq
# MF_3934_NoIndex_L002_R2_010.fastq
# MF_3934_NoIndex_L002_R2_020.fastq
# MF_3934_NoIndex_L002_R2_030.fastq
# MF_3934_NoIndex_L002_R2_040.fastq
# MF_3934_NoIndex_L002_R2_049.fastq
# MF_3934_NoIndex_L002_R2_050.fastq
# MF_3934_NoIndex_L002_R2_060.fastq
# MF_3934_NoIndex_L002_R2_061.fastq
# MF_4057_NoIndex_L006_R2_017.fastq
# MF_4057_NoIndex_L006_R2_043.fastq
# MF_4057_NoIndex_L006_R2_051.fastq
# MF_4057_NoIndex_L006_R2_052.fastq
# MF_4300_NoIndex_L004_R2_017.fastq
# MF_4300_NoIndex_L004_R2_033.fastq
# MF_4300_NoIndex_L004_R2_034.fastq
# MF_4300_NoIndex_L004_R2_041.fastq
# MF_4300_NoIndex_L004_R2_042.fastq
# MF_4300_NoIndex_L004_R2_050.fastq
# MF_4300_NoIndex_L004_R2_051.fastq
# MF_4306_NoIndex_L005_R2_009.fastq
# MF_4306_NoIndex_L005_R2_018.fastq
# MF_4306_NoIndex_L005_R2_027.fastq
# MF_4306_NoIndex_L005_R2_036.fastq
# MF_4306_NoIndex_L005_R2_045.fastq
# MF_4306_NoIndex_L005_R2_054.fastq
# MF_5269_NoIndex_L007_R2_032.fastq
# MF_5269_NoIndex_L007_R2_048.fastq

# 7  FAIL Per base sequence quality in D0LAC, all R1
# MF_4306_NoIndex_L005_R1_006.fastq
# MF_4306_NoIndex_L005_R1_007.fastq
# MF_4306_NoIndex_L005_R1_008.fastq
# MF_4306_NoIndex_L005_R1_017.fastq
# MF_4306_NoIndex_L005_R1_030.fastq
# MF_4306_NoIndex_L005_R1_032.fastq
# MF_4306_NoIndex_L005_R1_033.fastq

# (39+7)/1592 = 0.02889 
# However, things are not that serious. Go on align nevertheless. 

#######################################
# 2.2 novoalign use the illumina format fastq files

mkdir sam

for chip in D0LAC D0L9C; do 
    for sample in MF_1348 MF_2874 MF_3925 MF_3934 \
                  MF_4057 MF_4300 MF_4306 MF_5269; do
        qsub novo.sh $chip $sample 
    done
done

# because of using mpiexec, we now have all the mpiexec output in the sam files.
# use a python script to filt them. 

t1.py : 
#---------------------------------------
#!/home/klinbrc/bin/python3
import sys

ifile=open(sys.argv[1])
for line in ifile:
    if line.startswith('@') or line.startswith('HWI'):
        print(line[:-1])
ifile.close()
#---------------------------------------

for chip in D0LAC D0L9C; do 
    for sample in MF_1348 MF_2874 MF_3925 MF_3934 \
                  MF_4057 MF_4300 MF_4306 MF_5269; do
        ./t1.py sam/${chip}_${sample}.sam > t_${chip}_${sample}.sam  
        mv t_${chip}_${sample}.sam sam/${chip}_${sample}.sam
    done
done

# still got messages like "[proxy:0:3@node22] got pmi command (from 4): get"
# messed up in the sam file
# the mpi novoalign runs left these troubles. 

#######################################
# 2.3 convert to bam 

ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa

samtools view -bT $ref  sam/${chip}_${sample}.sam \
> sam/${chip}_${sample}.bam 

# last 12432, 11986 and 23867 lines are missing in 
# D0LAC_MF_4057.sam, D0L9C_MF_3925.sam and D0L9C_MF_4300.sam 
# because of the mpi messages...

#######################################
# 2.4 sort, remove duplications, indexing
samtools sort  sam/${chip}_${sample}.bam sam/${chip}_${sample}_sorted
samtools rmdup sam/${chip}_${sample}_sorted.bam sam/${chip}_${sample}_nodup.bam
samtools sort  sam/${chip}_${sample}_nodup.bam sam/${chip}_${sample}
samtools index sam/${chip}_${sample}.bam

#######################################
# 2.5 split into chr 
mkdir sam_chr 

for chip in D0LAC D0L9C; do 
    for sample in MF_1348 MF_2874 MF_3925 MF_3934 \
                  MF_4057 MF_4300 MF_4306 MF_5269; do
        for chr in 01 02 03 04 05 06 07 08 09 10 11 \
                   12 13 14 15 16 17 18 19 20 21 22 ; do 
            samtools view sam/${chip}_${sample}_nodup.bam chr${chr} -b \
                 > sam_chr/${chip}_${sample}_${chr}.bam 
            samtools index sam_chr/${chip}_${sample}_${chr}.bam 
        done
    done
done
# well, this is buggy. samltools chr will not recognize 01 02 etc. 

# with a s_chr = 1 2 3 etc for chr 01 02 03 etc 

~/bin/bamtools filter -in sam/${chip}_${sample}.bam -region ${s_chr} \
    -out sam_chr/${chip}_${sample}_${chr}.bam

#######################################
# 2.6 depth of coverage 

qsub depth.sh D0L9C MF_4306 

# or 
# download and unzip qualimap, change the line
# 

# JAVA_OPTS="-Xms32m -Xmx18G -XX:MaxPermSize=1024m"
# ~/bin/qualimap/qualimap bamqc -bam D0L9C_MF_1348_nodup.bam \
#  -outdir qualimap_D0L9C_MF_1348
# 
# it seems that qualimap asks for a lot of memory. 
# and it allways complains the input bam is not sorted according to the header.

# GATK works, 12 hours, 4G, gives us a summary only. 


################################################################################
# 3. What we got?                                                              #
################################################################################

8*2 samples, 100G fastq / sample, ~45G zipped.

~1000€ per sample
~500€ in two years.
# pipeline?

the ref has 3,101,804,739 sites
the fastq has 8,000,000 read (100 bp) * 42 ~ 61
the average coverage is 13 (12~14)
400000000
9,200,000

2679061*1.5
4,000,000 calls
3,500,000 SNPs

samtools view -F 4 sam/D0L9C_MF_1348_nodup.bam | perl -lane 'print "$F[2]\t$F[3]"' > t.hits
samtools view -F 4 sam/D0LAC_MF_5269_nodup.bam | perl -lane 'print "$F[2]\t$F[3]"' > t2.hits
cnv-seq.pl --test t.hits --ref ref.hits --genome  human > t.out
~/temp/cnv-seq/cnv-seq.pl --test t.hits --ref ref.hits --genome  human > t.out

################################################################################
# 4. CNV calling from BIC-seq                                                  #
################################################################################

# comparison between control and case 
# we have only cases so far 
# case1 vs case1 on different chips gives nothing 
# case1 vs case2 

install.packages('.', repos = NULL, type="source",lib="~/R")
library('BICseq', lib.loc="~/R")

bicseq <- BICseq(sample = "sam_chr/D0L9C_MF_5269_01.bam",
              reference = "sam_chr/D0LAC_MF_5269_01.bam", seqNames = c("1"));  
segs <- getBICseg(object = bicseq, bin = 100, lambda = 2, winSize = 200, quant = 0.95, mult = 1); 
seg.summary <- BICseq:::getSummary(segs,correction=TRUE);
write.table(seg.summary,'t.table',append=TRUE)
# etc...

bicseq <- BICseq(sample = "sam_chr/D0LAC_MF_5269_01.bam", 
              reference = "sam_chr/D0LAC_MF_4306_01.bam", seqNames = c("1")); 
segs <- getBICseg(object = bicseq, bin = 100, lambda = 2, winSize = 200, quant = 0.95, mult = 1); 
seg.summary <- BICseq:::getSummary(segs,correction=TRUE);
write.table(seg.summary,'t.table',append=TRUE)
# etc...


grep -v chr D0L9C_MF_5269_vs_D0LAC_MF_5269.table| grep -v Inf | \
 awk '{if ($2=="\"22\"") print $3,$7,"\n",$4,$7}' > t1.dat
grep -v chr D0LAC_MF_5269_vs_D0LAC_MF_4306.table| grep -v Inf | \
 awk '{if ($2=="\"22\"") print $3,$7,"\n",$4,$7}' > t2.dat 


################################################################################
# 5. CNVnator CNV calling                                                      #
################################################################################

# download root from http://root.cern.ch
# in ~/temp 
wget ftp://root.cern.ch/root/root_v5.34.04.source.tar.gz
tar xvzf root_v5.34.04.source.tar.gz

# pain to install 

# the bloody CNVnator v 0.2.7 doesn't work
# had to switch v 0.2.5
cd ~/temp
wget http://sv.gersteinlab.org/cnvnator/CNVnator_v0.2.5.zip
unzip CNVnator_v0.2.5.zip
cd CNVnator_v0.2.5/src
cd samtools 
make
cd ..
export LD_LIBRARY_PATH=~/local/root/lib:$LD_LIBRARY_PATH
export PATH=~/local/root/bin:$PATH
export ROOTSYS=~/local/root
make 

cd ~/dev/wgs

chr=22
sample=D0L9C_MF_5269

/home/kuang/temp/CNVnator_v0.2.5/src/cnvnator -root ${sample}_${chr}.root \
 -genome hg19 -chrom ${chr} -tree sam_chr/${sample}_${chr}.bam
 
# get the chr22.fa from /scratch/data/reference_genomes/human/human_g1k_v37.fa
# mv it to ref_seq dir

/home/kuang/temp/CNVnator_v0.2.5/src/cnvnator -root ${sample}_${chr}.root \
  -his 100 -d ref_seq

/home/kuang/temp/CNVnator_v0.2.5/src/cnvnator -root ${sample}_${chr}.root \
 -stat 100
/home/kuang/temp/CNVnator_v0.2.5/src/cnvnator -root ${sample}_${chr}.root \
 -partition 100
/home/kuang/temp/CNVnator_v0.2.5/src/cnvnator -root ${sample}_${chr}.root \
 -call 100 > ${sample}_${chr}.cnv

 
################################################################################
# 6. cnv-seq CNV calling                                                       #
################################################################################
samtools view -F 4 sam_chr/D0LAC_MF_5269_22.bam | perl -lane 'print "$F[2]\t$F[3]"' > D0LAC_MF_5269_22.hits
samtools view -F 4 sam_chr/D0L9C_MF_5269_22.bam | perl -lane 'print "$F[2]\t$F[3]"' > D0L9C_MF_5269_22.hits
samtools view -F 4 sam_chr/D0LAC_MF_4306_22.bam | perl -lane 'print "$F[2]\t$F[3]"' > D0LAC_MF_4306_22.hits

cnv-seq.pl --test D0LAC_MF_5269_22.hits --ref D0L9C_MF_5269_22.hits --genome-size 35194385
# give a file D0LAC_MF_5269_22.hits-vs-D0L9C_MF_5269_22.hits.log2-0.6.pvalue-0.001.minw-4.count
# had to download the file cnv.R
R
source('cnv.R')
ncnv<-cnv.cal("D0LAC_MF_5269_22.hits-vs-D0L9C_MF_5269_22.hits.log2-0.6.pvalue-0.001.minw-4.count", log2=0.6, min=4, chromosomal=TRUE, annota=TRUE);
write.table(ncnv,"ac_vs_9c.cnv",row.names=FALSE,sep="\t");

cnv-seq.pl --test D0LAC_MF_5269_22.hits --ref D0LAC_MF_4306_22.hits --genome-size 35194385
R

source('cnv.R')
ncnv<-cnv.cal("D0LAC_MF_5269_22.hits-vs-D0LAC_MF_4306_22.hits.log2-0.6.pvalue-0.001.minw-4.count", log2=0.6, min=4, chromosomal=TRUE, annota=TRUE);
write.table(ncnv,"52_vs_43.cnv",row.names=FALSE,sep="\t");

awk '{print $2,$3,$7,$8}' 52_vs_43.cnv | grep -v NA | grep -v log | \
    grep -v Inf | awk '{if ($4<0.00001)print $1,$3,"\n",$2,$3}' > t1.dat

awk '{print $2,$3,$7,$8}' ac_vs_9c.cnv | grep -v NA | grep -v log | \
    grep -v Inf | awk '{if ($4<0.00001)print $1,$3,"\n",$2,$3}' > t2.dat

# no CG content or whatever correction, it's not working properly

################################################################################
# 7. CNV from affy                                                             #
################################################################################
# WTCCCT535509 is MF_5269 col 220 
sort -k 1,1 \
    /home/kuang/dev/PE/wtccc2_pe_cel/CNV/GenomeWideSNP_6.hg18.cnv_defs > t.def
# mv the last line cnp_id etc to the first line 
paste t.def \
    /home/kuang/dev/PE/wtccc2_pe_cel/CNV/PE_affymetrix_CNV.canary_confs.txt \
    /home/kuang/dev/PE/wtccc2_pe_cel/CNV/PE_affymetrix_CNV.canary_calls.txt \
    > affy_.cnv

# 5819 x 2 + 1 + 1 + 4 = 11644
# 155 5975



################################################################################
#           CNV calling of the pilot data from Korean                          #
################################################################################

# working directory: /home/klinbrc/scratch/dev/cnv_pilot

# 13 individuals, 2 (_1/_2) fastq file per individual
# need to change ids: 
# mf_1391  to MF_1391 
# mf_1795  to MF_1795 
# mf_2864  to MF_2864
# mf_2874  to MF_2874 
# mf_3930  to MF_3930

# 13 x 2 fastq files in fastq/
 MF_1391 MF_1393 MF_1527 MF_1795 MF_2864 
 MF_2874 MF_3925 MF_3930 MF_3934 
 SGAP232 SGAP240 SGAP481 SRI1087 
# 2.17 T 

# fastqc is fine ...

# gzip them to save some space 
for ifile in fastq/* ; do 
    gzip $ifile 
done 

# or do it with a array job:

t.sh 
#---------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -cwd
#$-q long.q,bignode.q,short.q
#$-t 1-13
#$ -p -200

SAMPLE_LOOKUP=( MF_1391 MF_1393 MF_1527 MF_1795 MF_2864 MF_2874 \
        MF_3925 MF_3930 MF_3934 SGAP232 SGAP240 SGAP481 SRI1087 )

i=$(expr $SGE_TASK_ID - 1)

sample=${SAMPLE_LOOKUP[$i]}

gzip fastq/${sample}_1.fastq
gzip fastq/${sample}_2.fastq
#---------------------------------------
# end of array job t.sh 

qsub t.sh 

################################################################################
# 1. novoalign

# novocraft
# version 3.02.05 16 April 2014 
# I want the newest version because I can trust that the web manual should work.

# use reference genome /scratch/data/reference_genomes/human/human_g1k_v37.fa
# index it with novocraft

~/bin/novoindex hg19 /scratch/data/reference_genomes/human/human_g1k_v37.fa

# it uses at least 9G memory. 
# Done in 6 mins.

########################################
# 1.2 novoalign 

~/bin/novoalign -d hg19 -f fastq/MF_1391_1.fastq.gz fastq/MF_1391_2.fastq.gz

# novoalign uses about 10G memory. To be safe, allocate 12G for it. 
# The problem is we need to convert the sam into bam instantly to save space. 

novo.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -cwd
#$ -q long.q
#$ -l h_vmem=12G
#$ -pe multi_thread 4

sample=$1 

/share/apps/novocraft_20130909/bin/novoalign -c 4  \
    -d hg19 \
    -f fastq/${sample}_1.fastq.gz fastq/${sample}_2.fastq.gz \
    -o SAM $'@RG\tID:PSZ\tPL:ILLUMINA\tSM:'${sample} \
     | ~/bin/samtools view -S -b - | tee > ${sample}.bam
#--------------------------------------
# end of novo.sh 
     
for sample in  MF_1391 MF_1393 MF_1527 MF_1795 MF_2864 \
               MF_2874 MF_3925 MF_3930 MF_3934 \
               SGAP232 SGAP240 SGAP481 SRI1087 ; do 
    qsub novo.sh $sample 
done 

# two weeks per individual
# three weeks for MF_1527, maybe assigned to a slow node
# "multi_thread 4" and "-c 4 " there can be increased to 7 or maybe even 8

#######################################
# 1.3 sort, remove duplications, indexing

# backup the original bam files 
cp *.bam /archive/klinbrc/cnv_pilot_org_bam/

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -cwd
#$ -l h_vmem=6G
#$ -q long.q,bignode.q,short.q

sample=$1 

~/bin/samtools sort  ${sample}.bam           ${sample}_sorted
~/bin/samtools rmdup ${sample}_sorted.bam    ${sample}_nodup.bam

rm ${sample}.bam  ${sample}_sorted.bam

~/bin/samtools rmdup ${sample}_nodup.bam     ${sample}.bam 
~/bin/samtools index ${sample}.bam 
rm ${sample}_nodup.bam 
#--------------------------------------
# end of t.sh

for sample in  MF_1391 MF_1393 MF_1527 MF_1795 MF_2864 \
               MF_2874 MF_3925 MF_3930 MF_3934 \
               SGAP232 SGAP240 SGAP481 SRI1087 ; do 
    qsub t.sh $sample 
done 

########################################
# 1.4 split bam into chr bams

mkdir wg_bam
mv *.bam wg_bam/
mv *.bam.bai wg_bam/
mkdir chr_bam 

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -cwd
##$ -l h_vmem=6G
#$ -q long.q,bignode.q,short.q

sample=$1
chr=$2 

chr_int=$(echo $chr | sed 's/^[0]//')

~/bin/samtools view -bh wg_bam/$sample.bam $chr_int \
                     >  chr_bam/${sample}_${chr}.bam
~/bin/samtools index    chr_bam/${sample}_${chr}.bam

#--------------------------------------
# end of t.sh 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    for sample in  MF_1391 MF_1393 MF_1527 MF_1795 MF_2864 \
                   MF_2874 MF_3925 MF_3930 MF_3934 \
                   SGAP232 SGAP240 SGAP481 SRI1087 ; do 
        qsub t.sh $sample $chr 
    done 
done 

################################################################################
# 2 cnv freec callining                                                        #
################################################################################
cd /home/klinbrc/scratch/dev/wgs/cnv_freec
ls ../../cnv_pilot/chr_bam/*.bam > bam.ls

f_freec.py 
#--------------------------------------
#!/home/klinbrc/bin/python3

ifile=open('bam.ls')
bam_names=ifile.read().split()
ifile.close()

for bam_name in bam_names:
    chr=bam_name.split('_')[-1][:2]   # 21 
    conf_name='conf_'+bam_name.split('/')[4][:-4] #'conf_MF_1391_21'
    
    # write the conf file 
    conf_file=open(conf_name,'w')
    conf_file.write('[general]\nchrLenFile = hg19_chr'+chr+'.len\n')
    conf_file.write('window = 3000\nstep = 1000\nploidy = 2\n')
    conf_file.write('chrFiles = chr_fasta\n')
    conf_file.write('intercept=1\nminMappabilityPerWindow = 0.7\n')
    conf_file.write('outputDir = res \n\nsex=XX\nbreakPointType=4\n')
    conf_file.write('gemMappabilityFile = t2.mappability\n\n')
    conf_file.write('[sample]\n\n')
    conf_file.write('mateFile = '+bam_name+'\n\n')
    conf_file.write('inputFormat = BAM\nmateOrientation = FR\n\n')
    conf_file.write('[control]\n\n[BAF]\n')
    conf_file.close()
    
#--------------------------------------
# end of run_freec.py 

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -cwd
#$ -l h_vmem=3G
#$ -q long.q,bignode.q,short.q

PATH=.:~/bin:$PATH
export PATH

sample=$1
chr=$2 

~/bin/freec -conf conf_${sample}_$chr
#--------------------------------------
# end of t.sh 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do
    for sample in  MF_1391 MF_1393 MF_1527 MF_1795 MF_2864 \
                   MF_2874 MF_3925 MF_3930 MF_3934 \
                   SGAP232 SGAP240 SGAP481 SRI1087 ; do 
        qsub t.sh $sample $chr 
    done 
done 

################################################################################
# 3 NimbleGen CNV calling                                                      #
################################################################################

# availabe at cgh/
