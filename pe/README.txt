################################################################################
#                        QC and GWAS of the PE data                            #
################################################################################

1. the original data. 

 evokerSet
 
 The directory of the intensity files. 
 Evoker can use this directory to show the genotyping calling errors.
 
 evokerSet/PE_affymetrix.sample
 
 The original set has 5831 individuals. 196 individuals have one Sanger ids but
 two different (plate) Ids. The plate with more genotype calling will be used.

 
################################################################################
# 2. Initial QC                                                                #
################################################################################

 2.1 The genotype data files were converted to plink format using GTOOL 
 with the default threshold of 0.9 on p values.
 
 pe.bed
 pe.fam
 pe.bim
 
 5635 individuals, 929556 SNPs.
 
################################################################################
 2.2 sex check
 
 pe1.bed
 pe1.fam
 pe1.bim
 
 58 individuals failed plink sex check were removed.
 
################################################################################
 2.3 remove X, Y, XY, MT chrs
 
 pe2.bed
 pe2.fam
 pe2.bim
 
 38895 SNPs removed.
 5577 individuals, 890661 SNPs left
 
################################################################################
# 3. family structure                                                          #
################################################################################
 
 fam/fam.zsh

################################################################################
# 4. QC                                                                        #
################################################################################
 
 qc/qc.zsh
 
################################################################################
# 5. PCA                                                                       #
################################################################################
 
 pca/pca.zsh
 
################################################################################
# 6. subset                                                                    #
################################################################################

 subset/subset.zsh
 
################################################################################
# 7. address the meeting note, on pca vs meta analysis                         #
################################################################################

 meta/meta.zsh

################################################################################
# 8. alternative phenotype                                                     #
################################################################################

 pheno/pheno.zsh
 
