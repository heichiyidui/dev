
################################################################################
# rat differential expression analysis from sequencing data                    #
################################################################################

################################################################################
# input:
# rat cortical neuron cells 
# illumina sequencing data
# 4 abeta, 4 dkk1, 3 wildtype control (1 control failed QC)

# top differentially expressed genes? 
# alternative splicing in the gene? 

################################################################################
# 1. fastqc                                                                    #
################################################################################

# working directory
cd /home/klinbrc/scratch/dev/rat_rna

for ifile in data/Sample_RK*/*.gz; do gunzip  $ifile; done
for ifile in data/Sample_RK*/*.gz; do qsub gu.sh $ifile; done
for ifile in data/Sample_RK*/*.fastq; do fastqc $ifile; done 

rm -r */*/*_fastqc
mv */*/*_fastqc.zip .

# all failed Per base sequence content (That's normal. It happens a lot.)
# all failed Per base GC content	(That happens. Don't worry.)

# The per base sequence quality is pretty good. 

################################################################################

#### the protocol ########
wget http://www.bioconductor.org/help/course-materials/2013/CSAMA2013/\
tuesday/afternoon/DESeq_protocol.pdf

###########################
# bioconductor install 
sudo R 
source("http://bioconductor.org/biocLite.R")
biocLite()
install.packages('RColorBrewer')
biocLite(c('ggplot2','scales'))

###########################
# install libxml2 and gcc-c++ in yast 
biocLite(c("annotate"))
biocLite(c("DESeq", "edgeR"))
biocLite(c("ShortRead", "GenomicRanges"))
biocLite("org.Dm.eg.db")

###########################
# install HTSeq 
# install python-numpy python-numpy-devel python-matplotlib via yast
wget https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.5.4p5.tar.gz
tar xvzf HTSeq-0.5.4p5.tar.gz 
cd HTSeq-0.5.4p5
python setup.py build
sudo python setup.py install

###########################
# tophat2 install 
wget http://tophat.cbcb.umd.edu/downloads/tophat-2.0.10.Linux_x86_64.tar.gz
tar xvzf tophat-2.0.10.Linux_x86_64.tar.gz 
mv tophat-2.0.10.Linux_x86_64 ~/bin
cd ~/bin
ln -s tophat-2.0.10.Linux_x86_64/tophat2 .
ln -s tophat-2.0.10.Linux_x86_64/tophat .

cd /home/klinbrc/scratch/dev/rat_rna
rm tophat-2.0.10.Linux_x86_64.tar.gz

###########################
#bowtie2 install 
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/\
bowtie2-2.1.0-linux-x86_64.zip/download
mv download bowtie2-2.1.0-linux-x86_64.zip    
unzip bowtie2-2.1.0-linux-x86_64.zip
mv bowtie2-2.1.0 ~/bin
rm bowtie2-2.1.0-linux-x86_64.zip
cd ~/bin 
ln -s bowtie2-2.1.0/bowtie2 .
ln -s bowtie2-2.1.0/bowtie2-align .
ln -s bowtie2-2.1.0/bowtie2-build .
ln -s bowtie2-2.1.0/bowtie2-inspect .

cd /home/klinbrc/scratch/dev/rat_rna

################################################################################
# tophat 2 align 

# download the rat ref (the cluster ref is not permitted!)
mkdir ref 
cd ref
wget http://hgdownload.soe.ucsc.edu/goldenPath/rn5/bigZips/rn5.fa.gz
# (version 09-May-2012)
gunzip rn5.fa.gz
# build bowtie index 
bowtie2-build rn5.fa rn5
/share/bin/bowtie-build rn5.fa rn5

# gene annotation information
# wget ftp://ftp.ensembl.org/pub/release-74/gtf/rattus_norvegicus/\
# Rattus_norvegicus.Rnor_5.0.74.gtf.gz
# the gtf file I got from ftp.ensembl.org is of wrong format. 
# go http://genome.ucsc.edu/cgi-bin/hgTables?command=start
# track Ensembl genes, use ensGene table. 
# choose GTF as ouput format. 
# (Thanks to Damian Kao)



# run tophat

sloc=data/Sample_RK01-C1/RK01-C1_ATCACG
id=$(echo $sloc| awk -F"-" '{print substr($2,1,2)}')

gfile=ref/rn5.gtf

tophat -G $gfile -o $id ref/rn5  \
${sloc}_L006_R1_001.fastq,${sloc}_L007_R1_001.fastq \
${sloc}_L006_R2_001.fastq,${sloc}_L007_R2_001.fastq 

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m e
#$ -N tophat2
#$ -cwd
#$ -q long.q,bignode.q
#$ -l h_vmem=7G

sloc=$1
# sample location like "data/Sample_RK01-C1/RK01-C1_ATCACG"
id=$(echo $sloc| awk -F"-" '{print substr($2,1,2)}')
# sample id like "C1"

gfile=ref/rn5.gtf

PATH=~/bin:$PATH

~/bin/tophat -p 4 -G $gfile -o $id ref/rn5  \
${sloc}_L006_R1_001.fastq,${sloc}_L007_R1_001.fastq \
${sloc}_L006_R2_001.fastq,${sloc}_L007_R2_001.fastq 

#--------------------------------------

for sloc in data/Sample_RK01-C1/RK01-C1_ATCACG \
            data/Sample_RK02-C2/RK02-C2_CGATGT \
            data/Sample_RK03-C3/RK03-C3_TTAGGC \
            data/Sample_RK05-A1/RK05-A1_ACAGTG \
            data/Sample_RK06-A2/RK06-A2_GCCAAT \
            data/Sample_RK07-A3/RK07-A3_CAGATC \
            data/Sample_RK08-A4/RK08-A4_ACTTGA \
            data/Sample_RK09-D1/RK09-D1_GATCAG \
            data/Sample_RK10-D2/RK10-D2_TAGCTT \
            data/Sample_RK11-D3/RK11-D3_GGCTAC \
            data/Sample_RK12-D4/RK12-D4_CTTGTA ; do 
    qsub t.sh $sloc 
done 
# takes forever
# DO NOT USE THE SHORT QUEUE!

# with -p 4 for the multi-thread running (4 times faster)
# the jobs were killed on the cluster before the last step 
# had to create the last accepted_hits.bam via another script

t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m e
#$ -N smt 
#$ -cwd
#$ -q long.q,bignode.q
#$ -l h_vmem=4G
PATH=~/bin:$PATH 

i_id=$1
samtools sort ${i_id}/tmp/accepted_hits0.bam ${i_id}/tmp/accepted_hits0_sorted
samtools sort ${i_id}/tmp/accepted_hits1.bam ${i_id}/tmp/accepted_hits1_sorted
samtools sort ${i_id}/tmp/accepted_hits2.bam ${i_id}/tmp/accepted_hits2_sorted
samtools sort ${i_id}/tmp/accepted_hits3.bam ${i_id}/tmp/accepted_hits3_sorted

samtools merge -h ${i_id}/tmp/rn5_genome.bwt.samheader.sam \
    ${i_id}/accepted_hits.bam ${i_id}/tmp/accepted_hits?_sorted.bam 

/home/klinbrc/bin/tophat-2.0.10.Linux_x86_64/bam_merge -Q \
    --sam-header ${i_id}/tmp/rn5_genome.bwt.samheader.sam \
    ${i_id}/unmapped.bam \
    ${i_id}/tmp/unmapped_left_0.bam ${i_id}/tmp/unmapped_right_0.bam

#--------------------------------------

################################################################################
# sorting (via gene name or coordinates) and indexing via samtools
t.sh 
#--------------------------------------
#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m e
#$ -N t_sort
#$ -cwd
#$ -q long.q,bignode.q,short.q

source ~/.bashrc

i_id=$1
samtools sort -n ${i_id}/accepted_hits.bam ${i_id}_sn
samtools sort    ${i_id}/accepted_hits.bam ${i_id}_s

samtools view -o ${i_id}_sn.sam ${i_id}_sn.bam
samtools index ${i_id}_s.bam
#--------------------------------------

t.sh A4 
# 1.1G 1.5 hours

################################################################################
# get the reading depths
i_id=A4 
for i_id in A1 A2 A3 A4 C1 C2 C3 D1 D2 D3 D4 ; do 
    htseq-count -s no -a 10 ${i_id}_sn.sam ref/rn5.gtf > ${i_id}.count
done 
# around 4 hours each

################################################################################
# EdgeR simple analysis                                                        #
################################################################################

sample.csv
#--------------------------------------
ID,condition,LibraryLayout
A1,ABE,double
A2,ABE,double
A3,ABE,double
A4,ABE,double
C1,CTL,double
C2,CTL,double
C3,CTL,double
D1,DKK,double
D2,DKK,double
D3,DKK,double
D4,DKK,double
#--------------------------------------

# in R
library("edgeR")

samples=read.csv("sample.csv", stringsAsFactors=FALSE)
samples$countf = paste(samples$ID, "count", sep=".")

###########################
# read counts and QC
counts = readDGE(samples$countf)$counts

noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual", 
           "not_aligned","alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms>1)>=3 & !noint
# sum(counts) is 690511233
# nrow(counts) is 29193
counts = counts[keep,]
# sum(counts) is 78830537
# nrow(counts) is 10816

############################
# grouping 
d = DGEList(counts=counts, group=samples$condition)
d = calcNormFactors(d)

plotMDS(d, labels=samples$shortname, 
    col=c("darkgreen","blue")[factor(samples$condition)])
# the PCA like plot is not that promising

d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
# looks fine
plotBCV(d)
# looks fine


de1 = exactTest(d, pair=c('ABE','CTL')) # comparison between ABE and CTL
de2 = exactTest(d, pair=c('CTL','DKK')) # comparison between CTL and DKK

tt1 = topTags(de1, n=nrow(d)) # nothing is significant according to FDR
tt2 = topTags(de2, n=nrow(d)) # even worse!

# tt1 top ENSRNOT00000024752 is from CCSER1 
# Diseases associated with CCSER1 include schizophrenia.

# tt2 top ENSRNOT00000043827 is from ST8SIA3
# Diseases associated with ST8SIA3 include glioblastoma, and among its related 
# super-pathways are WNT ligand biogenesis and trafficking and N-Glycan antennae
# elongation.

################################
# DESeq

samples=read.csv("sample.csv", stringsAsFactors=FALSE)
samples$countf = paste(samples$ID, "count", sep=".")

samplesDESeq = with(samples, data.frame(
    shortname = I(ID),
    countf = I(countf),
    condition = condition, 
    LibraryLayout = LibraryLayout))

library("DESeq")
cds = newCountDataSetFromHTSeqCount(samplesDESeq)

# normalization
cds = estimateSizeFactors(cds)
sizeFactors(cds)

# pca
cdsB = estimateDispersions(cds, method="blind")
vsd = varianceStabilizingTransformation(cdsB)
plotPCA(vsd, intgroup="condition")
# still the PCA plot is messy 

cds = estimateDispersions(cds)
plotDispEsts(cds)
# this plot is fine

res1 = nbinomTest(cds,"CTL","ABE")
res2 = nbinomTest(cds,"CTL","DKK")
plotMA(res1)
plotMA(res2)

# the P values are even larger

# res1 top ENSRNOT00000024752 is the same 
# res2 top ENSRNOT00000043227 is not identified. 
# ENSRNOT00000043827 is #26 in res2 

################################################################################
# last bits 

# get the gene ids from ensMart 
http://www.ensembl.org/Multi/martview
# get the corresponding networks from genemania 
http://www.genemania.org/

# A-beta vs ctrl has some hits.
# dkk vs ctrl got nothing. 

################################################################################
# remove the outliers and re-run the edge analysis                             #
################################################################################

# remove 4 and 11 

sample.csv
#--------------------------------------
ID,condition,LibraryLayout
A1,ABE,double
A2,ABE,double
A3,ABE,double
C1,CTL,double
C2,CTL,double
C3,CTL,double
D1,DKK,double
D2,DKK,double
D3,DKK,double
#--------------------------------------

# in R
library("edgeR")

samples=read.csv("sample.csv", stringsAsFactors=FALSE)
samples$countf = paste(samples$ID, "count", sep=".")

###########################
# read counts and QC
counts = readDGE(samples$countf)$counts

noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual", 
           "not_aligned","alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms>1)>=3 & !noint
# sum(counts) is 563262773
# nrow(counts) is 29193
counts = counts[keep,]
# sum(counts) is 61095722
# nrow(counts) is 10589

############################
# grouping 
d = DGEList(counts=counts, group=samples$condition)
d = calcNormFactors(d)

plotMDS(d, labels=samples$shortname, 
    col=c("darkgreen","blue")[factor(samples$condition)])
# the PCA like plot slightly better.

d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
# looks fine
plotBCV(d)
# looks fine

de1 = exactTest(d, pair=c('ABE','CTL')) # comparison between ABE and CTL
de2 = exactTest(d, pair=c('CTL','DKK')) # comparison between CTL and DKK

tt1 = topTags(de1, n=nrow(d)) # nothing is significant according to FDR
tt2 = topTags(de2, n=nrow(d)) # even worse!

# tt1 top ENSRNOT00000066191 is from Runx1t1 
# FDR is way better! Even significant!

# tt2 top ENSRNOT00000075366 is from a pseudogene
# FDR is better. Still not significant. 

################################################################################

