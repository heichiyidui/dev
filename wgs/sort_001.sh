#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m beas
#$ -N sort1
#use current directory as working directory
#$ -cwd

#$ -q long.q,bignode.q
# Declare how much memory is required PER slot - default is 2Gbytes
#$ -l h_vmem=6G
################################################################################

sample=D0LAC_MF_1348_001_uns.bam
describer=$(echo ${sample} | sed 's/.bam//' | sed 's/_uns//') 

samdir=/home/klinbrc/scratch/dev/wgs/sam
picarddir=/share/apps/picard-tools_1.35/jar

/usr/java/latest/bin/java -Xmx4g -jar ${picarddir}/SortSam.jar \
SO=coordinate VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \
TMP_DIR=${samdir}/tmp \
MAX_RECORDS_IN_RAM=250000 \
INPUT=${samdir}/${sample} \
OUTPUT=${samdir}/${describer}_sorted.bam \


/usr/java/latest/bin/java -Xmx4g -jar ${picarddir}/MergeSamFiles.jar \
SO=coordinate VALIDATION_STRINGENCY=SILENT \
TMP_DIR=tmp \

RGID=String
RGPL=illumina \
RGPU=
ASSUME_SORTED=true \

OUTPUT=File \
INPUT=File \