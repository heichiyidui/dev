#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m beas
#$ -N samcall
#use current directory as working directory
#$ -cwd

#$ -q long.q,bignode.q
# Declare how much memory is required PER slot - default is 2Gbytes
#$ -l h_vmem=10G
################################################################################

ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa 

samtools mpileup -uf $ref  D0LAC_MF_1348_*_nodup.bam \
 | bcftools view -bvcg - > D0LAC_MF_1348.raw.bcf  

bcftools view D0LAC_MF_1348.raw.bcf \
 | vcfutils.pl varFilter -D 100 > D0LAC_MF_1348.flt.vcf 
