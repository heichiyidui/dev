#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
#$ -q bignode.q,long.q
#$ -l h_vmem=8G
#$ -o /scratch/home/ylimbrc/haloplex/cluster_msg/4_sam2bam/output
#$ -e /scratch/home/ylimbrc/haloplex/cluster_msg/4_sam2bam/error


Sample_Barcode=$1

infile=/scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/${Sample_Barcode}_L00?.sam
filepath=$(echo $infile | sed 's/.sam//' )

java \
 -Xmx6g \
 -Djava.io.tmpdir=/tmp \
 -jar ~/bin/picard-tools-1.113/SortSam.jar \
 SO=coordinate \
 INPUT=${filepath}.sam \
 OUTPUT=${filepath}.bam \
 VALIDATION_STRINGENCY=LENIENT \
 CREATE_INDEX=true
