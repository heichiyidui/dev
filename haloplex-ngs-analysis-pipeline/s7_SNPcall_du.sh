#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
#$ -q bignode.q,long.q
#$ -l h_vmem=8G
#$ -o /scratch/home/ylimbrc/haloplex/cluster_msg/7_SNPcall/output
#$ -e /scratch/home/ylimbrc/haloplex/cluster_msg/7_SNPcall/error


Sample_Barcode=$1

infile=/scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/${Sample_Barcode}_L00?.marked.realigned.bam
filepath=$(echo $infile | sed 's/.marked.realigned.bam//' )

/share/java/jdk1.7.0/bin/java \
 -Xmx6g \
 -jar ~/bin/GenomeAnalysisTK.jar \
 -T UnifiedGenotyper \
 -R /scratch/home/ylimbrc/haloplex/refgenome/hg19.fa \
 --dbsnp /scratch/home/ylimbrc/haloplex/dbsnp/00-All.vcf \
 -stand_call_conf 50.0 \
 -stand_emit_conf 10.0 \
 -dcov 1000 \
 -I ${filepath}.marked.realigned.bam \
 -o ${filepath}.snp.raw.vcf

