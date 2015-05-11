#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
#$ -q bignode.q,long.q
#$ -l h_vmem=8G
#$ -o /scratch/home/ylimbrc/haloplex/cluster_msg/1_trimming/output
#$ -e /scratch/home/ylimbrc/haloplex/cluster_msg/1_trimming/error

Sample_Barcode=$1

Sample=$(echo $Sample_Barcode | sed 's/_[ATGC]*//' )
adapseq=$(echo $Sample_Barcode | sed 's/DCR[0-9]*_//' )
mkdir /scratch/home/ylimbrc/haloplex/Trimmed/${Sample_Barcode}

~/bin/trim_galore \
    --phred33 \
    -o /scratch/home/ylimbrc/haloplex/Trimmed/${Sample_Barcode} \
    --paired \
    --clip_R1 5 \
    --clip_R2 5 \
    --adapter ${adapseq} \
    /scratch/home/ylimbrc/haloplex/P335/demultiplexed/Project_Hodges_335/Sample_${Sample}/${Sample_Barcode}_L00?_R1_001.fastq.gz \
    /scratch/home/ylimbrc/haloplex/P335/demultiplexed/Project_Hodges_335/Sample_${Sample}/${Sample_Barcode}_L00?_R2_001.fastq.gz

