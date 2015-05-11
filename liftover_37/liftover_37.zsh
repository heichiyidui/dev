
################################################################################
# To liftover SNPs chr and positions to build 37                               #
# To flip all SNPs to the + strand                                             #
################################################################################

################################################################################
# 1. the original plink files 

tar xvf LIFTOVER_37.tar
rm LIFTOVER_37.tar

################################################################################
# 2. I hate long file names 

mkdir origin_bims

cp LIFTOVER_37/GWAS1raw_Ammar.bim                                      origin_bims/g1.bim 
cp LIFTOVER_37/GWAS2raw_Ammar.bim                                      origin_bims/g2.bim
cp LIFTOVER_37/irish.bim                                               origin_bims/ir.bim 
cp LIFTOVER_37/M_Depr_UK_controls_cleaned_SNPs_samples_Ammar_plink.bim origin_bims/uc.bim
cp LIFTOVER_37/meta3_noNIH.bim                                         origin_bims/mt.bim 
cp LIFTOVER_37/MNDFull_alleleTop.bim                                   origin_bims/mn.bim
cp LIFTOVER_37/NIH_Corriel.bim                                         origin_bims/co.bim 
cp LIFTOVER_37/NIH_IT_ACGT.bim                                         origin_bims/it.bim
cp LIFTOVER_37/pheALS.bim                                              origin_bims/ph.bim

################################################################################
# 3. the strand files 

mkdir strand_files
cd strand_files

# from the website http://www.well.ox.ac.uk/~wrayner/strand/
#660
wget http://www.well.ox.ac.uk/~wrayner/strand/Human660W-Quad_v1_A-b37.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/Human660W-Quad_v1_C-b37.zip
#610
wget http://www.well.ox.ac.uk/~wrayner/strand/Human610Quadv1_B-b37-strand-v2.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/Human610-Quadv1_B-b37-SourceStrand.zip
#370
wget http://www.well.ox.ac.uk/~wrayner/strand/HumanCNV370v1_C-b37-strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/HumanCNV370v1_C-b37-SourceStrand.zip
#550
wget http://www.well.ox.ac.uk/~wrayner/strand/BDCHP-1X10-HUMANHAP550_11218540_C-b37-strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/HumanHap550-2v3_B-b37-strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/HumanHap550-2v3_B-b37-SourceStrand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/HumanHap550-2v3_B-b37-Illmn.strand.zip
#300
wget http://www.well.ox.ac.uk/~wrayner/strand/BDCHP-1x10-HUMANHAP300v1-1_11219278_C-b37-strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/HumanHap300v2_A-b37-strand.zip

for zfile in *.zip; do unzip -o $zfile ; done 
# some .miss and .multiple files have the same names, we overwrite them with -o

# we got 12 strand files 

# BDCHP-1x10-HUMANHAP300v1-1_11219278_C-b37.strand
# BDCHP-1X10-HUMANHAP550_11218540_C-b37.strand
# Human610-Quadv1_B-b37-SourceStrand.strand
# Human610-Quadv1_B-b37-v2.strand
# Human660W-Quad_v1_A-b37.strand
# Human660W-Quad_v1_C-b37.strand
# HumanCNV370v1_C-b37-SourceStrand.strand
# HumanCNV370v1_C.strand
# HumanHap300v2_A-b37.strand
# HumanHap550-2v3_B-b37.Illmn.strand
# HumanHap550-2v3_B-b37-SourceStrand.strand
# HumanHap550-2v3_B-b37.strand

cd ..

################################################################################
# 4. check the overlapping of SNP ids between .bim and .strand files

# a python script will do 
chk_id_overlap.py

# run it like 
chk_id_overlap.py origin_bims/g1.bim strand_files/HumanHap550-2v3_B-b37.strand

# the output will be 
317503 561329 6183 311320 250009
# from left to right 
# the number of SNP ids in the bim file
# the number of SNP ids in the strand file 
# the number of SNP ids in the bim but not the strand file
# the number of overlapping SNP ids (in both files)
# the number of SNP ids in the strand but not the bim file 

# use for loops to run it over 9 bim file vs 12 strand files 

rm t.out
for bfile in origin_bims/*.bim ; do 
    for sfile in strand_files/*.strand; do
        chk_id_overlap.py $bfile $sfile >> t.out
    done
done

# to check the order of the bim files 
for bfile in origin_bims/*.bim ; do 
    echo $bfile
done

# to check the order of the strand files 
for sfile in strand_files/*.strand; do
    echo $sfile
done

# see liftover_37.xls

################################################################################
# 5. check the allele of overlapping SNPs

# The strand files have the illumina top strand SNP allele
# check the SNPs in the bim file have the same allele 

# a python script chk_allele.py
# run it like 
chk_allele.py origin_bims/g1.bim strand_files/HumanHap550-2v3_B-b37.strand

# the output will be 
311320 155834

# the number of overlapping SNPs and 
# in these SNPs the number of SNPs with the same allele

# In the strand file Human610-Quadv1_B-b37-SourceStrand.strand
# we don't have the allele!

mv strand_files/Human610-Quadv1_B-b37-SourceStrand.strand \
   strand_files/Human610-Quadv1_B-b37-SourceStrand.strand.bak


# use for loops to run it over 9 bim file vs 11 strand files 

rm t.out
for bfile in origin_bims/*.bim ; do 
    for sfile in strand_files/*.strand; do
        chk_allele.py $bfile $sfile >> t.out
    done
done






  


