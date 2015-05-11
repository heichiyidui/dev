#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m beas
#$ -N qualimap1
#use current directory as working directory
#$ -cwd

#$ -q long.q,bignode.q
# Declare how much memory is required PER slot - default is 2Gbytes
#$ -l h_vmem=13G
################################################################################

chip=D0L9C
sample=MF_1348

JAVA_OPTS="-Xms32m -Xmx12G -XX:MaxPermSize=1024m"
~/bin/qualimap/qualimap bamqc -bam ${chip}_${sample}_nodup.bam \
  -outdir qualimap_${chip}_${sample}
