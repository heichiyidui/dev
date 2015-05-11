#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
#$ -q bignode.q,long.q
#$ -l h_vmem=8G
#$ -o /scratch/home/ylimbrc/haloplex/cluster_msg/3_sai2sam/output
#$ -e /scratch/home/ylimbrc/haloplex/cluster_msg/3_sai2sam/error


Sample_Barcode=$1

infile=/scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/R1/${Sample_Barcode}_L00?_R1_001.sai
outfile=$(echo $infile | sed 's/sai/sam/' | sed 's/_R1_001//' | sed 's/\/R1//' )
#filename=$(echo $infile | sed 's/\/scratch\/home\/ylimbrc\/haloplex\/Aligned\/'${Sample}'\/R1\///' | sed 's/_R1_001//' | sed 's/.sai//' )
#id=$(echo $filename | sed 's/_L00[1-8]*//' )
sm=$(echo $Sample_Barcode | sed 's/_[AGTC]*//' )

~/bin/bwa sampe \
 -r "@RG\tID:${Sample_Barcode}\tLB:hg19\tSM:${sm}\tPL:ILLUMINA" \
 /scratch/home/ylimbrc/haloplex/refgenome/hg19 \
 /scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/R1/${Sample_Barcode}_L00?_R1_001.sai \
 /scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/R2/${Sample_Barcode}_L00?_R2_001.sai \
 /scratch/home/ylimbrc/haloplex/Trimmed/${Sample_Barcode}/${Sample_Barcode}_L00?_R1_001_val_1.fq.gz \
 /scratch/home/ylimbrc/haloplex/Trimmed/${Sample_Barcode}/${Sample_Barcode}_L00?_R2_001_val_2.fq.gz \
 > $outfile
