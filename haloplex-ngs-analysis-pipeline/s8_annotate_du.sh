#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
#$ -q bignode.q,long.q
#$ -l h_vmem=8G
#$ -o /scratch/home/ylimbrc/haloplex/cluster_msg/8_annotate/output
#$ -e /scratch/home/ylimbrc/haloplex/cluster_msg/8_annotate/error


Sample_Barcode=$1

infile=/scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/${Sample_Barcode}_L00?.snp.raw.vcf
filepath=$(echo $infile | sed 's/.snp.raw.vcf//' )

~/bin/annovar/convert2annovar.pl \
 --format vcf4 \
 --includeInfo \
 ${filepath}.snp.raw.vcf > ${filepath}.snp.annovar

~/bin/annovar/summarize_annovar.pl \
 --buildver hg19 \
 ${filepath}.snp.annovar \
 /scratch/home/ylimbrc/haloplex/humandb \
 -outfile ${filepath}.snps
