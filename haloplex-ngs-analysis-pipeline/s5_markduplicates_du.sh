#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
#$ -q bignode.q,long.q
#$ -l h_vmem=8G
#$ -o /scratch/home/ylimbrc/haloplex/cluster_msg/5_markduplicates/output
#$ -e /scratch/home/ylimbrc/haloplex/cluster_msg/5_markduplicates/error


Sample_Barcode=$1

infile=/scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/${Sample_Barcode}_L00?.bam
filepath=$(echo $infile | sed 's/.bam//' )

java \
 -Xmx6g \
 -Djava.io.tmpdir=/tmp \
 -jar ~/bin/picard-tools-1.113/MarkDuplicates.jar \
 INPUT=${filepath}.bam \
 OUTPUT=${filepath}.marked.bam \
 METRICS_FILE=metrics \
 CREATE_INDEX=true \
 VALIDATION_STRINGENCY=LENIENT
