#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
#$ -q bignode.q,long.q
#$ -l h_vmem=8G
#$ -o /scratch/home/ylimbrc/haloplex/cluster_msg/6_indelrealign/output
#$ -e /scratch/home/ylimbrc/haloplex/cluster_msg/6_indelrealign/error


Sample_Barcode=$1

infile=/scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/${Sample_Barcode}_L00?.marked.bam
filepath=$(echo $infile | sed 's/.marked.bam//' )

/share/java/jdk1.7.0/bin/java \
 -Xmx6g \
 -jar ~/bin/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -R /scratch/home/ylimbrc/haloplex/refgenome/hg19.fa \
 -I ${filepath}.marked.bam \
 -o ${filepath}.bam.list

/share/java/jdk1.7.0/bin/java \
 -Xmx6g \
 -Djava.io.tmpdir=/tmp \
 -jar ~/bin/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R /scratch/home/ylimbrc/haloplex/refgenome/hg19.fa \
 -targetIntervals ${filepath}.bam.list \
 -I ${filepath}.marked.bam \
 -o ${filepath}.marked.realigned.bam
