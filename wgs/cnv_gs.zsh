################################################################################
# CNV calling via GENOME STRIP                                                 #
################################################################################

################################################################################
# INPUT                                                                        #
################################################################################

cd /home/klinbrc/scratch/dev/cnv_gs/

# 1. Analysis ready BAM files 
#  - Aligned, sorted, indexed, duplicates removed.
#  - from the wgs.zsh script 

ls sam_chr/
# 8x2 individuals, 22 chrs, 352 bam files and 352 bam index files

# need to add the library tag into the bam files 
# otherwise genome strip complains about ids like:
# Error: Cannot determine library identifier for read HWI-ST1168:61:D0LACAC ... 
# should go back to wgs.zsh for to add the lib tag during alignment

# see sam_tag.sh 
# for bamfile in sam_chr/*.bam ; do qsub sam_tag.sh $bamfile ; done

# 2. reference
export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa

# 3. Alignability Mask
# - Indicates which reference positions are uniquely alignable
# - Must be based on the same reference you are using
# - Commonly used masks are available for download
# from ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/
# human_g1k_v37.mask.100.fasta  1000 Genomes Phase I (build 37) read length 100

mkdir align_mask
cd align_mask
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/human_g1k_v37.mask.100.fasta.gz 
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/human_g1k_v37.mask.100.fasta.sizes
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/human_g1k_v37.mask.100.fasta.fai
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/human_g1k_v37.mask.100.fasta.stats
  gunzip human_g1k_v37.mask.100.fasta.gz
cd ..
# small files

# 4 CN2 mask
# - Flags regions unlikely to be copy number polymorphic, used for estimate GC bias
# - CN2 masks for common reference sequences are available for download
# from ftp://ftp.broadinstitute.org/pub/svtoolkit/cn2masks/

mkdir cn2_mask
cd cn2_mask
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/cn2masks/cn2_mask_g1k_v37.fasta.fai
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/cn2masks/cn2_mask_g1k_v37.fasta.gz
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/cn2masks/cn2_mask_g1k_v37.mask.bed
  gunzip cn2_mask_g1k_v37.fasta.gz
cd ..
# small files 

# genderMapFile etc are ignored. We are not dealing with sex chrs yet. 

################################################################################
# PREPROCESSING                                                                #
################################################################################

mkdir tmp      # temp directory 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    mkdir metadata_${chr}     # the output directories of preprocessing 
    ls sam_chr/*_${chr}.bam > bam_${chr}.list 
done 

export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export SV_DIR=/share/apps/genomestrip_1.03
# SV_DIR should be a system variable after the installation of Genome STRiP
export classpath=${SV_DIR}/lib/SVToolkit.jar:\
${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
java -Xmx4g -cp ${classpath}  \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVPreprocess.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${classpath} \
    -configFile ${SV_DIR}/conf/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -R $ref \
    -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
    -md metadata_${chr} \
    -I bam_${chr}.list \
    -run 
done 

# rm -r metadata_??/* tmp/* Q*bignode*.out

# This is running on bignode alone.
# around 7 hours running on chr 01 

# when used with qsub 
# many jobs got error message with 
# Exit code: 134 
# 
# guess it can be done with '-T' etc 
# but not with '-S ... ?????.q'

################################################################################
# SVDiscovery                                                                  #
################################################################################

export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export SV_DIR=/share/apps/genomestrip_1.03

export classpath=${SV_DIR}/lib/SVToolkit.jar:\
${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
java -Xmx4g -cp ${classpath} \
    org.broadinstitute.sv.main.SVDiscovery \
    -T SVDiscovery \
    -configFile ${SV_DIR}/conf/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -R $ref \
    -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
    -md metadata_${chr} \
    -I bam_${chr}.list \
    -O ${chr}.sites.vcf \
    -runDirectory metadata_${chr} 
done 

# 0.5 hour on chr 22, about 27 hours on chr 01
# could save some running time via setting maximumSize

# no additional filters used 

# this can 
################################################################################
# alt allele alignment                                                         #
################################################################################

export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export SV_DIR=/share/apps/genomestrip_1.03

export classpath=${SV_DIR}/lib/SVToolkit.jar:\
${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
java  -Xmx4g -cp ${classpath} \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVAltAlign.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -configFile ${SV_DIR}/conf/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -cp ${classpath} \
    -R $ref \
    -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
    -md metadata_${chr} \
    -runDirectory metadata_${chr} \
    -vcf ${chr}.sites.vcf \
    -I bam_${chr}.list \
    -O ${chr}.deletions.alt.bam \
    -run
done 

# eventually, no bam produced because of ??.sites.alt.fasta being empty files 
# no breakpoint found? 

################################################################################
# SV genotyping                                                                #
################################################################################
export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export SV_DIR=/share/apps/genomestrip_1.03

export classpath=${SV_DIR}/lib/SVToolkit.jar:\
${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
java -Xmx4g -cp ${classpath} \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVGenotyper.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${classpath} \
    -configFile ${SV_DIR}/conf/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -R $ref \
    -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
    -runDirectory metadata_${chr} \
    -md metadata_${chr} \
    -I bam_${chr}.list \
    -vcf ${chr}.sites.vcf \
    -O ${chr}.genotypes.vcf \
    -run 
done 

################################################################################
# comparison between two plates                                                #
################################################################################

#######################################
# 1. preprocessing 

mkdir tmp      # temp directory 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    for plate in D0LAC D0L9C ; do 
        mkdir metadata_${plate}_${chr} # the output directories of preprocessing
        ls sam_chr/${plate}_*_${chr}.bam > ${plate}_${chr}_bam.list 
    done 
done 

export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export SV_DIR=/share/apps/genomestrip_1.03
export classpath=${SV_DIR}/lib/SVToolkit.jar:\
${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    for plate in D0LAC D0L9C ; do 
java -Xmx4g -cp ${classpath}  \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVPreprocess.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${classpath} \
    -configFile ${SV_DIR}/conf/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -R $ref \
    -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
    -md metadata_${plate}_${chr} \
    -I  ${plate}_${chr}_bam.list  \
    -run 
    done 
done 
# finished in 2 days 

#######################################
# 2. SV Discovery     

mkdir res 

export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export SV_DIR=/share/apps/genomestrip_1.03
export classpath=${SV_DIR}/lib/SVToolkit.jar:\
${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    for plate in D0LAC D0L9C ; do 
    
java -Xmx4g -cp ${classpath} \
    org.broadinstitute.sv.main.SVDiscovery \
    -T SVDiscovery \
    -configFile ${SV_DIR}/conf/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -R $ref \
    -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
    -md metadata_${plate}_${chr} \
    -I  ${plate}_${chr}_bam.list \
    -O res/${plate}_${chr}.sites.vcf \
    -runDirectory metadata_${plate}_${chr}
     
    done 
done 

#######################################
# 3. SV genotyping       

export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export SV_DIR=/share/apps/genomestrip_1.03
export classpath=${SV_DIR}/lib/SVToolkit.jar:\
${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    for plate in D0LAC D0L9C ; do 
java -Xmx4g -cp ${classpath} \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVGenotyper.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${classpath} \
    -configFile ${SV_DIR}/conf/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -R $ref \
    -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
    -runDirectory metadata_${plate}_${chr} \
    -md metadata_${plate}_${chr} \
    -I  ${plate}_${chr}_bam.list \
    -vcf res/${plate}_${chr}.sites.vcf  \
    -O res/${plate}_${chr}.genotypes.vcf \
    -run 
    done 
done 
################################################################################


