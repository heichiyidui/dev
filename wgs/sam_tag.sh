#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m e
#$ -N sam_tag
#use current directory as working directory
#$ -cwd

#$ -q long.q,bignode.q,short.q
#------------------------------------------------------------------------------#

bamfile=$1
samfile=${bamfile/.bam/.sam}

samtools view -h $bamfile > $samfile

sed -i 's/ID:PSZ\tPL:ILLUMINA/ID:PSZ\tPL:ILLUMINA\tLB:tlib/' $samfile
sed -i 's/PG:Z:novoalignMPI\tRG:Z:PSZ/PG:Z:novoalignMPI\tRG:Z:PSZ\tLB:Z:tlib/' \
    $samfile

samtools view -bS $samfile > $bamfile 
samtools index    $bamfile 

rm $samfile 
