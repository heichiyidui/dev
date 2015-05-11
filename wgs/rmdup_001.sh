#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m beas
#$ -N rmdup1
#use current directory as working directory
#$ -cwd

#$ -q long.q,bignode.q
# Declare how much memory is required PER slot - default is 2Gbytes
#$ -l h_vmem=8G
################################################################################

sample=D0LAC_MF_1348_001_sorted.bam
describer=$(echo ${sample} | sed 's/.bam//' | sed 's/_sorted//') 

samdir=/home/klinbrc/scratch/dev/wgs/sam
picarddir=/share/apps/picard-tools_1.35/jar

/usr/java/latest/bin/java -Xmx6g -jar ${picarddir}/MarkDuplicates.jar \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
TMP_DIR=${samdir}/tmp \
INPUT=${samdir}/${sample} \
OUTPUT=${samdir}/${describer}_nodup.bam \
METRICS_FILE=${samdir}/${describer}_metric.file
