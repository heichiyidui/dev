################################################################################
# CNV calling with PINDEL                                                      #
################################################################################

################################################################################
# install 
cd ~/temp
git clone git://github.com/genome/pindel.git
cd pindel
./INSTALL ~/bin/samtools-0.1.19/
# the samtools source code dir can also be
# /share/apps/samtools_0.1.18/samtools-0.1.18/

cd demo 
time ../pindel -f simulated_reference.fa -i simulated_config.txt -o output
# vs 
time pindel -f simulated_reference.fa -i simulated_config.txt -o output

# This Pindel (version 0.2.5a1, July 23 2013) is faster than 
# the system default /share/bin/pindel (version 0.2.5, June 4 2013)
# user    23.500s vs user    41.010s

# I hope the memory leaking problem was solved...

# Pindel version 0.2.5a8, Oct. 14 2014.
cd ~/temp 
mv pindel ~/bin/pindel_0.2.5a8 

cd ~/bin 
ln -s ~/bin/pindel_0.2.5a8/pindel .
 
# install the sam2pindel tool for sam files produced by aligner other than BWA
cd /home/klinbrc/bin/pindel_0.2.5a8/src

g++ -O3 SAM_2_PINDEL_cin_2011Sept26.cpp -o sam2pindel
cd ~/bin
ln -s ~/bin/pindel_0.2.5a8/src/sam2pindel .

################################################################################
export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa

for sample in MF_2874 MF_3925 MF_3934; do 
    samtools view /isilon/panfs/data/KuangBam/1312KHS-0010/${sample}_sorted.bam\
       | sam2pindel - Output4Pindel.txt 300 $sample 0 Illumina-PairEnd
    pindel -f $ref -p Output4Pindel.txt -c ALL -o $sampe 
done 

export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa
export bam_dir=/isilon/panfs/data/KuangBam/1312KHS-0010

for sample in MF_2874 MF_3925 MF_3934; do 
    for chr in 01 02 03 04 05 06 07 08 09 10 11 \
               12 13 14 15 16 17 18 19 20 21 22 ; do 
        chr_int=$(echo $chr | sed 's/^[0]//')
    
        samtools view $bam_dir/${sample}_$chr.bam\
         | sam2pindel - pin.txt 400 $sample 0 Illumina-PairEnd
        sed -i "s/chr$chr_int/$chr_int/" pin.txt 
        pindel -f $ref -p pin.txt -c $chr_int -o ${sample}_$chr 
    done
done



# prepare the configuration files
for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    ls sam_chr/*_${chr}.bam > t.ls 
    awk -F"_" '{print $3 "_" $4}' t.ls > t2.ls
    paste t.ls t2.ls | awk '{print $1,400,$2}' > ${chr}.cfg
done 
rm t.ls t2.ls 

# use 400 as the expected insert size 

pindel -f $ref -i 01.cfg -c  1 -o o_01
pindel -f $ref -i 02.cfg -c  2 -o o_02
pindel -f $ref -i 03.cfg -c  3 -o o_03
pindel -f $ref -i 04.cfg -c  4 -o o_04
pindel -f $ref -i 05.cfg -c  5 -o o_05
pindel -f $ref -i 06.cfg -c  6 -o o_06
pindel -f $ref -i 07.cfg -c  7 -o o_07
pindel -f $ref -i 08.cfg -c  8 -o o_08
pindel -f $ref -i 09.cfg -c  9 -o o_09
pindel -f $ref -i 10.cfg -c 10 -o o_10
pindel -f $ref -i 11.cfg -c 11 -o o_11
pindel -f $ref -i 12.cfg -c 12 -o o_12
pindel -f $ref -i 13.cfg -c 13 -o o_13
pindel -f $ref -i 14.cfg -c 14 -o o_14
pindel -f $ref -i 15.cfg -c 15 -o o_15
pindel -f $ref -i 16.cfg -c 16 -o o_16
pindel -f $ref -i 17.cfg -c 17 -o o_17
pindel -f $ref -i 18.cfg -c 18 -o o_18
pindel -f $ref -i 19.cfg -c 19 -o o_19
pindel -f $ref -i 20.cfg -c 20 -o o_20
pindel -f $ref -i 21.cfg -c 21 -o o_21
pindel -f $ref -i 22.cfg -c 22 -o o_22

# 4G to 71G memory, toooooo much 

################################################################################
# do it on each individual

#$ -tc 40

export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    chr_int=$(echo $chr | sed 's/^[0]//')
    for sample in 1348 2874 3925 3934 4057 4300 4306 5269; do 
        for plate in D0LAC D0L9C ; do 
            ls sam_chr/${plate}_MF_${sample}_${chr}.bam | \
                awk -F"_" '{print $0,400,"MF_" $4}' > t.cfg 
            pindel -f $ref -i t.cfg -c $chr_int -o ${plate}_MF_${sample}_${chr}
        done 
    done 
done 
# 3.5g memory on chr 22, 3.7g on chr 01 

#for chr in 01 02 03 04 05 06 07 08 09 10 11 \
#           12 13 14 15 16 17 18 19 20 21 22 ; do 
for chr in  01  ; do 
    chr_int=$(echo $chr | sed 's/^[0]//')
    for sample in 1348 2874 3925 3934 4057 4300 4306 5269; do 
        for plate in D0L9C D0LAC; do 
             qsub p.sh $chr $chr_int $sample $plate 
        done 
    done 
done 
# so far no problem. 20G memory is fine for chr 12 - 22
# using 30G for 07-11 
# using 40G for 04-06
# using 50G for 01-03 

# often got error messages like 
# Loading reference genome ...
# Error: fasta line starts with ^@ instead of '>'. Aborting.

# resubmit the job (many times) might slove it. 
# don't submit those jobs via loop 




export ref=/scratch/data/reference_genomes/human/human_g1k_v37.fa

for chr in 01 02 03 04 05 06 07 08 09 10 11 \
           12 13 14 15 16 17 18 19 20 21 22 ; do 
    chr_int=$(echo $chr | sed 's/^[0]//')
    for sample in MF_2874 MF_3925 MF_3934; do 
        echo "/isilon/panfs/data/KuangBam/1312KHS-0010/"\
${sample}_${chr}.bam 400 $sample > t.cfg 
        #pindel -f $ref -i t.cfg -c $chr_int -o ${sample}_${chr}
    done 
done 

~/bin/pindel -f $ref -i t.cfg -c ALL -o $sample
