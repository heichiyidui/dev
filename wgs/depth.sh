#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m beas
#$ -N depth1
#use current directory as working directory
#$ -cwd

#$ -q long.q,bignode.q
# Declare how much memory is required PER slot - default is 2Gbytes
#$ -l h_vmem=6G
################################################################################

chip=$1
sample=$2

ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa 
samdir=/home/klinbrc/scratch/dev/wgs/sam
gatk=/share/apps/GenomeAnalysisTK_1.0.3471/GenomeAnalysisTK.jar

samtools index ${chip}_${sample}_nodup.bam 

/usr/java/latest/bin/java -Xmx3g -jar $gatk \
-T DepthOfCoverage \
-R $ref \
--omitDepthOutputAtEachBase \
--omitLocusTable \
-I ${samdir}/${chip}_${sample}_nodup.bam  \
-o ${samdir}/${chip}_${sample}.DepthOfCoverage \
