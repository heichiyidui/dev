################################################################################
# CNV calling via DELLY                                                        #
################################################################################

################################################################################
# 1. get the program from http://www.embl.de/~rausch/delly.html

cd ~/temp 
http://www.embl.de/~rausch/delly_v0.0.11.tar.gz

git clone --recursive https://github.com/tobiasrausch/delly.git
# something this connection is not working on the brc cluster

# wget https://codeload.github.com/tobiasrausch/delly/tar.gz/v0.6.1

tar xvzf delly_v0.0.11.tar.gz
cd delly_v0.0.11/
mv duppy delly invy jumpy ~/bin 

################################################################################
# 2

/share/bin/delly -o del.vcf \
   -g /scratch/data/reference_genomes/human/human_g1k_v37.fa \
   /isilon/panfs/data/KuangBam/1312KHS-0010/MF_3925_22.bam 
