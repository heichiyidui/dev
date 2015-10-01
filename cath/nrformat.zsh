rm nr.??.tar.gz

wget -r -nd -A "nr.??.tar.gz" ftp://ftp.ncbi.nlm.nih.gov/blast/db

# 14 files at 17/09/2013
# 15 files at 22/01/2014
# 17 files at 03/03/2014
# 17 files at 11/05/2014
# 26 files at 06/10/2014
# 28 files at 10/01/2015 
# 40 files at 30/09/2015 

for nfile in nr.??.tar.gz; do 
    tar xvzf $nfile
done

rm nr.??.tar.gz

#blastdbcmd -db nr -entry all -out nr.fasta -outfmt %f
#
#pfilt nr.fasta > nrfilt.fasta
#
#rm nr.fasta
#
#makeblastdb -title nrfilt -in nrfilt.fasta -out nrfilt -dbtype prot
#
#rm nrfilt.fasta
