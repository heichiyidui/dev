################################################################################
# CNV calling via FREEC (Control-FREEC)                                        #
################################################################################

################################################################################
# 1. installation                                                              #
################################################################################

#######################################
# 1.1 the executable 
mkdir freec_src
cd freec_src/
  wget http://bioinfo-out.curie.fr/projects/freec/src/FREEC_Linux64.tar.gz
  tar xvzf FREEC_Linux64.tar.gz
  make 
  # had some warnings but works fine with gcc version 4.1.2 20080704
  # failed to compile with gcc 4.7.2 20130108
  # ThreadPool.h:101:22: error: 'usleep' was not declared in this scope
  # add '#include <unistd.h>' into ThreadPool.h and it compiles
  mv freec ~/bin  
  cd ..
rm -r freec_src/

#######################################
# 1.2 scripts 
mkdir freec_scripts
cd freec_scripts
  wget http://bioinfo-out.curie.fr/projects/freec/src/scripts_FREEC.zip
  unzip scripts_FREEC.zip 
  cd ..

#######################################
# 1.3 the chromosomes length files 
export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
perl freec_scripts/get_fasta_lengths.pl $ref 

head -n 22 res_human_g1k_v37.fa  | \
 awk '{ print $1 "\tchr" NR "\t" $4}' | split -l 1 

mv xaa hg19_chr01.len
mv xab hg19_chr02.len
mv xac hg19_chr03.len
# ...
mv xau hg19_chr21.len
mv xav hg19_chr22.len

mkdir chr_len
mv hg19_chr??.len chr_len

#######################################
# 1.4 chrFile chromosomes 
# to split the reference fasta file into 22 chromosomes
mkdir chr_fasta

# in python3
#--------------------------------------
#!/usr/bin/env python3
import fastaio
seqs=fastaio.read_fasta('/scratch/data/reference_genomes/human/human_g1k_v37.fa')
for i in range(22):
    fastaio.write_fasta('chr_fasta/chr'+str(i+1)+'.fasta',seqs[i:i+1])
#--------------------------------------

# use up to 16G memory, don't try this anywhere other than bignode 

#######################################
# 1.5 the mappablity file 

# download GEM from http://algorithms.cnag.cat/wiki/The_GEM_library, then
bunzip2 GEM-binaries-Linux-x86_64-core_2-20130406-045632.tbz
tar xvf GEM-binaries-Linux-x86_64-core_2-20130406-045632.tar
cd GEM-binaries-Linux-x86_64-core_2-20130406-045632/bin
# complains about the kernel being too old on the brc cluster. 
# had to use a desktop PC for it. 

# to get the index 
gem-indexer -i $ref -o t1 
# 7G memory, ~5 hours

# to get the Mappability file 
gem-mappability  -I t1.gem -l 100 -o t2 
# 7G memory, ~40 hours 
cd ../../
mv GEM-binaries-Linux-x86_64-core_2-20130406-045632/bin/t2.mappability .

################################################################################
# 2. run                                                                       #
################################################################################

mkdir res 

ls sam_chr/*.bam > bam.ls

run_freec.py 
#--------------------------------------
#!/usr/bin/env python3

ifile=open('bam.ls')
bam_names=ifile.read().split()
ifile.close()

for bam_name in bam_names:
    
    chr=bam_name.split('_')[-1][:2]               #'09'
    conf_name='conf_'+bam_name[:10]               #'MF_1234_22'
    
    # write the conf file 
    conf_file=open(conf_name,'w')
    conf_file.write('[general]\n')
    conf_file.write('chrLenFile = chr_len/hg19_chr'+chr+'.len\n')
    conf_file.write('window = 3000\nstep = 1000\nploidy = 2\n')
    conf_file.write('chrFiles = chr_fasta\n')
    conf_file.write('intercept=1\nminMappabilityPerWindow = 0.7\n')
    conf_file.write('outputDir = res \n\nsex=XX\nbreakPointType=4\n')
    conf_file.write('gemMappabilityFile = t2.mappability\n\n')
    conf_file.write('[sample]\n\n')
    conf_file.write('mateFile = ')
    conf_file.write('/isilon/panfs/data/KuangBam/1312KHS-0010/'+bam_name+'\n\n')
    conf_file.write('inputFormat = BAM\nmateOrientation = FR\n\n')
    conf_file.write('[control]\n\n[BAF]\n')
    conf_file.close()
    
    from subprocess import call
    call('~/bin/freec -conf ./'+conf_name,shell=True)
    call('rm ./'+conf_name,shell=True)
    #break;

#--------------------------------------

# 3 samples in 5 hours 


