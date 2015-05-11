#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m beas
#$ -N sam
#use current directory as working directory
#$ -cwd

#$ -q long.q,bignode.q
# Declare how much memory is required PER slot - default is 2Gbytes
#$ -l h_vmem=20G
################################################################################

ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa 
samdir=/home/klinbrc/scratch/dev/wgs/sam
picarddir=/share/apps/picard-tools_1.35/jar

for sample in *_novo.sam ; do 
    describer=$(echo ${sample} | sed 's/.sam//' | sed 's/_novo//') 
    echo $describer
    
    # Convert file from SAM to BAM format  
    samtools view -bT $ref $sample > ${describer}_uns.bam
    
#    rm $sample
done

# merge
samtools merge D0L9C_MF_1348_all.bam *_uns.bam

rm *_uns.bam

# sort 
/usr/java/latest/bin/java -Xmx12g -jar ${picarddir}/SortSam.jar \
SO=coordinate VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \
TMP_DIR=${samdir}/tmp \
MAX_RECORDS_IN_RAM=250000 \
INPUT=${samdir}/D0L9C_MF_1348_all.bam \
OUTPUT=${samdir}/D0L9C_MF_1348_sorted.bam \

# 2.2g memory for two chunks and lots of temp files!

rm D0L9C_MF_1348_all.bam

# remove duplicates
/usr/java/latest/bin/java -Xmx12g -jar ${picarddir}/MarkDuplicates.jar \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
TMP_DIR=${samdir}/tmp \
INPUT=${samdir}/D0L9C_MF_1348_sorted.bam \
OUTPUT=${samdir}/D0L9C_MF_1348_nodup.bam \
METRICS_FILE=${samdir}/D0L9C_MF_1348_metric.file

rm D0L9C_MF_1348_sorted.bam

# Depth of coverage

/usr/java/latest/bin/java -Xmx12g \
-jar /share/apps/GenomeAnalysisTK_1.0.3471/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $ref \
-I ${samdir}/D0L9C_MF_1348_nodup.bam \
-o ${samdir}/D0L9C_MF_1348.DepthOfCoverage \
--outputFormat csv -ct 1 -ct 2


