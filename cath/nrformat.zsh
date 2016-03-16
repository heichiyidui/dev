################################################################################
# NCBI nr database

rm nr.??.tar.gz

wget -r -nd -A "nr.??.tar.gz" ftp://ftp.ncbi.nlm.nih.gov/blast/db

# 14 files at 17/09/2013
# 15 files at 22/01/2014
# 17 files at 03/03/2014
# 17 files at 11/05/2014
# 26 files at 06/10/2014
# 28 files at 10/01/2015
# 40 files at 30/09/2015
# 48 files at 16/03/2016

for nfile in nr.??.tar.gz; do
    tar xvzf $nfile
done

rm nr.??.tar.gz

#######################################
# filter nr database?

# Remove various non-globular/biased regions using David Jone's pfilt program

# blastdbcmd -db nr -entry all -out nr.fasta -outfmt %f
#
# pfilt nr.fasta > nrfilt.fasta
#
# rm nr.fasta
#
# makeblastdb -title nrfilt -in nrfilt.fasta -out nrfilt -dbtype prot
#
# rm nrfilt.fasta

################################################################################
# UniRef 90 database

# maybe uniref90 is more appropriate.
# A. It's smaller. Jobs are faster.
# B. It has taxonomy. Multiple alignment mapping for prediction of protein
# interactions might be much easier.
# Downloaded uniref90 17/02/16
wget \
ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz

# 9.2G file, not too big at all.

gunzip uniref90.fasta.gz

~/bin/ncbi-blast-2.3.0+/bin/makeblastdb \
    -title uniref90 \
    -in uniref90.fasta \
    -out uniref90 \
    -dbtype prot

# 40253516 sequences
# about 48% of the 83476574 sequences of the nr database
################################################################################
