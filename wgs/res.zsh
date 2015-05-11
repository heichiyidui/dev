################################################################################
#   show the results                                                           #
################################################################################

cd res 

################################################################################
# 1. sample 
# MF_1348 at FlowCell_D0L9C

# WTCCCT535088 the 255th col in affy files 

################################################################################
# 2. affy 6.0 calling. 

# 1316 CNVs
# 7657804 callings, 3.3% 0, 8.3% 1, 85% 2, 3.6% 3 and more 

# use the QC threshold 0.01 as mentioned in the README.txt
# removed 621 of 1316 calling 

# remove 19 on X and Y chromosomals

# sort -g -k 1,2 -s 

# 676 CNV callings 

#     27 0
#     56 1
#    581 2
#      6 3
#      5 4
#      1 5

# 95 CNVs 
# most are deletions

# saved into in affy_.cnv 

# converting to build 37, using UCSC liftover, failed 9
##Partially deleted in new #chr1	146780429	147056085
##Split in new             #chr2	110173879	110519409
##Deleted in new           #chr7	141969866	142005168
##Split in new             #chr7	141693868	141712586
##Split in new             #chr7	141891589	141902885
##Split in new             #chr7	142155609	142167486
##Partially deleted in new #chr9	38906782	65412427

# CNV sizes: min 669, max 0.5M, mean 25k, median 8k

awk '{if ($1==1) print $2-1,2,"\n",$2,$4,"\n",$3,$4,"\n",$3+1,2}' affy_.cnv |\
    sort -g > affy.dat 
grace affy.dat 

################################################################################
# 3. Pindel 

# small in/dels based on break points 
# supposted to be very fast suffix tree algorithm 
# break point indentification might benefit from multiple samples 
# but too many samples slow it down a lot!
# chr 18 still running after 13 days!

# pindel finds are small 
# at least two reads support 

# small insertions (SI)
awk '{print $3 }' pindel/o_??_SI > t.in
listdis t.in > t.out 
# size:      650345
# min:       1.0000
# max:       82.0000
# mean:      4.5472
# std:       7.7942
# quartiles: 1 1 2 4 82.0000
# 8 tiles:   1 1 1 1 2 3 4 8 82.0000


awk '{if ($3>47) print $10-1,2,"\n",$10,4,"\n",$11,4,"\n",$11+1,2}' \
    pindel/o_01_SI > t2.dat
grace affy.dat t2.dat

# large insertions (LI)
# nothing found there

# deletions 
awk '{print $3 }' pindel/o_??_D > t.in
listdis t.in > t.out 
# size:      845722
# min:       1.0000
# max:       153069334.0000
# mean:      3777.9145
# std:       389557.5245
# quartiles: 1 1 2 4 153069334
# 8 tiles:   1 1 1 1 2 2 4 10 153069334

awk '{print $16}' pindel/o_01_D | h # number of supporting reads 

awk '{if ($3>22500) print $10-1,2,"\n",$10,0,"\n",$11,0,"\n",$11+1,2}' \
    pindel/o_01_D > t2.dat
grace affy.dat t2.dat


################################################################################
# 4. freec 

# sliding window, read counts, GC correction, 
# very very very fast ...
# output is smiple enough 

# insertions 
awk '{if ($4>2) print $3-$2}' freec/* > t.in
listdis t.in > t.out 

# size:      337
# min:       2999 
# max:       21205999 
# mean:      98521 
# std:       1152635 
# quartiles: 2999 6999 16999 46999 21205999 
# 8 tiles:   2999 5999 6999 10999 16999 29999 46999 76999 21205999 

# deletions

awk '{if ($4<2) print $3-$2}' freec/* > t.in
listdis t.in > t.out 

# size:      637
# min:       4999 
# max:       18377999 
# mean:      125610 
# std:       751672 
# quartiles: 4999 8999 35999 98999 18377999.0000
# 8 tiles:   4999 7999 8999 19999 35999 57999 98999 203999 18377999 

awk '{print $2-1,2,"\n",$2,$4,"\n",$3,$4,"\n",$3+1,2}' \
    freec/D0L9C_MF_1348_01.bam_CNVs > t2.dat 
grace affy.dat t2.dat &


################################################################################
# 5. genomestrip 

# deletions only. insertions promissed but never implemented. 

# 4007 deletions
grep -v # MF_1348_MF_1348_MF_1348_*.vcf | awk '{print $8}' |  head -n 1

# length GSDOPT
grep -v # gs_geno/MF_1348_MF_1348_MF_1348_*.vcf | awk '{print $8}' |  \
  awk -F";" '{print $17}' | awk -F"=" '{print $2}' > t.in 

listdis t.in > t.dat
# size:      4007
# min:       0.0000
# max:       239713767
# mean:      9527574
# std:       29192989.0391
# quartiles: 0 2121 11821 388172 239713767
# 8 tiles:   0 1072 2121 4635 11821 39687 388172 12832036 239713767

  
grep -v # gs_geno/MF_1348_MF_1348_MF_1348_01.genotypes.vcf | awk '{print $8}' |\
    awk -F";" '{print $3,$17}' | sed 's/\=/\ /g' | \
    awk  '{print $2-$4-1,2,"\n",$2-$4,0,"\n",$2,0,"\n",$2+1,2}' > t2.dat 

grace affy.dat t2.dat 

grep -v # gs_geno/MF_1348_MF_1348_MF_1348_01.genotypes.vcf | awk '{print $8}' | 

    



