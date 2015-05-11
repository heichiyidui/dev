#!/bin/sh
#$ -S /bin/sh
#$ -cwd 
#$ -q bignode.q,long.q
#$ -l h_vmem=8G
#$ -o /scratch/home/ylimbrc/haloplex/cluster_msg/2_alignment/output
#$ -e /scratch/home/ylimbrc/haloplex/cluster_msg/2_alignment/error


Sample_Barcode=$1

for Read in R1 R2 ; do
    #Sample=$(echo $Sample_Barcode | sed 's/_[ATGC]*//' )
    mkdir -p /scratch/home/ylimbrc/haloplex/Aligned/${Sample_Barcode}/${Read}
    infile=/scratch/home/ylimbrc/haloplex/Trimmed/${Sample_Barcode}/${Sample_Barcode}_L00?_${Read}_001_val_?.fq.gz
    outfile=$(echo $infile | sed 's/fq.gz/sai/' | sed 's/_val_[1-2]*//' | sed 's/Trimmed/Aligned/' | sed 's/\/'${Sample_Barcode}'\//\/'${Sample_Barcode}'\/'${Read}'\//' )

    ~/bin/bwa aln \
     -t 4 \
     -f ${outfile} \
     /scratch/home/ylimbrc/haloplex/refgenome/hg19 \
     /scratch/home/ylimbrc/haloplex/Trimmed/${Sample_Barcode}/${Sample_Barcode}_L00?_${Read}_001_val_?.fq.gz
done
