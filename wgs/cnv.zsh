################################################################################
# CNV calling                                                                  #
################################################################################

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

sam_chr/
# 8x2 individuals, 22 chrs, 352 bam files and 352 bam index files

# 2. reference
export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa

# 3. Alignability Mask
# - Indicates which reference positions are uniquely alignable
# - Must be based on the same reference you are using
# - Commonly used masks are available for download
# from ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/
# human_g1k_v37.mask.100.fasta  1000 Genomes Phase I (build 37) read length 100

  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/human_g1k_v37.mask.100.fasta.gz 
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/human_g1k_v37.mask.100.fasta.sizes
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/human_g1k_v37.mask.100.fasta.fai
  wget ftp://ftp.broadinstitute.org/pub/svtoolkit/svmasks/human_g1k_v37.mask.100.fasta.stats
mkdir align_mask
mv human_g1k_v37.mask.100.fasta.* align_mask
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

cd ..

gunzip */*.fasta.gz

################################################################################
# PREPROCESSING                                                                #
################################################################################

# need to add the library tag into the bam files 
# otherwise genome strip complains about ids like:
# Error: Cannot determine library identifier for read HWI-ST1168:61:D0LACAC ... 

# should go back to wgs.zsh for to add the lib tag during alignment

# see sam_tag.sh 

for bamfile in sam_chr/*.bam ; do qsub sam_tag.sh $bamfile ; done

# now, run the preprocessing set 
export gs_home=/share/apps/genomestrip_1.03
export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa

mkdir tmp      # temp directory 

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    mkdir metadata_${chr}     # the output directory for preprocessing 
    ls sam_chr/*_${chr}.bam > bam_${chr}.list 
done 

# rm -r metadata_??/* tmp/* Q*bignode*.out

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    java -Xmx4g -cp ${gs_home}/lib/gatk/Queue.jar:\
${gs_home}/lib/gatk/GenomeAnalysisTK.jar:\
${gs_home}/lib/SVToolkit.jar\
     org.broadinstitute.sting.queue.QCommandLine \
    -S ${gs_home}/qscript/SVPreprocess.q \
    -S ${gs_home}/qscript/SVQScript.q \
    -md metadata_${chr} \
    -gatk ${gs_home}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${gs_home}/lib/gatk/Queue.jar:\
${gs_home}/lib/gatk/GenomeAnalysisTK.jar:\
${gs_home}/lib/SVToolkit.jar\
    -tempDir tmp \
    -R $ref \
    -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
    -I bam_${chr}.list \
    -run &
done 

# This is running on bignode alone. Could spread the jobs to other nodes. 
# around 7 hours running on chr 01 

################################################################################
# SVDiscovery                                                                  #
################################################################################

export gs_home=/share/apps/genomestrip_1.03
export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export SV_DIR=/share/apps/genomestrip_1.03
# the variable SV_DIR should be a system variable after the installation of 
# GenomeSTRiP

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    java -Xmx4g -cp \
${gs_home}/lib/gatk/GenomeAnalysisTK.jar:\
${gs_home}/lib/SVToolkit.jar\
        org.broadinstitute.sv.main.SVDiscovery \
        -T SVDiscovery \
        -configFile ${gs_home}/conf/conf/genstrip_parameters.txt \
        -md metadata_${chr} \
        -R $ref \
        -genomeMaskFile align_mask/human_g1k_v37.mask.100.fasta \
        -I bam_${chr}.list \
        -O ${chr}.sites.vcf \
        -runDirectory metadata_${chr}
done

# 0.5 hour on chr 22, more than 27 hours on chr 01

################################################################################
# FREEC (Control-FREEC)                                                        #
################################################################################

# 1. installation is very easy 
wget http://bioinfo-out.curie.fr/projects/freec/src/FREEC_Linux64.tar.gz
tar xvzf FREEC_Linux64.tar.gz
make 
mv freec ~/bin  # version v6.4 

# 2. download relavent files 
# hs19_chr19.len
# GC_profile.cnp
# hg19_snp131.SingleDiNucl.1based.txt
wget http://xfer.curie.fr/get/7hZIk1C63h0/hg19_len100bp.tar.gz

# 3. run 
freec -conf config_chr19.txt 
# very fast

################################################################################
# PINDEL                                                                       #
################################################################################

/share/bin/pindel
# Pindel version 0.2.5, June 4 2013
pindel -f $ref \
 -i D0LAC_MF_5269_19_pindel.cfg \
 -o D0LAC_MF_5269_19 \
 -c 19 
# 4G to 12G memory 


