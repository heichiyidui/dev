################################################################################
#                      CKB batch effect (manual) check                         #
################################################################################

################################################################################
# Given batches of SNP calling from AxiomGT1, regressions have be performed to #
# detect SNPs with potential batch and other effects.                          #
# Some manual (visual) checking might be neccessary for deciding thresholds    #
# and manual removal of SNPs with bad calling.                                 #
################################################################################

cd /kuser/kuangl/dev/ckb_batch_check

################################################################################
# 1. Input files

#######################################
# 1.1 SNP lists:

# on the nc2 server:
# plate effect:
/kuser/shared/data/GWAS_backup/full_data/plate-effect/\
variant_plate_effects_v2.txt
# 33621 entries, 30570 uniq SNPs

# batch effect:
/kuser/shared/data/GWAS_backup/full_data/batch_test/\
variant_batch_effects.txt
# 6407 entries, 4048 uniq SNPs

#######################################
# 1.2 calling files:

# on the nc2 server again:
# 7 batches of calling files at
/kuser/shared/data/GWAS_backup

plates1-53
plates54-105
plates106-156
plates157-209
plates210-261
plates262-318
plates319-367

# From each batch, we need four files:
AxiomGT1.calls.txt
AxiomGT1.confidences.txt
AxiomGT1.snp-posteriors.txt
AxiomGT1.summary.txt

################################################################################
# 2. SNP cluster plots

mkdir b01 b02 b03 b04 b05 b06 b07

tail -n +2 /kuser/shared/data/GWAS_backup/full_data/batch_test/\
variant_batch_effects.txt | awk '{print $2}' | sort | uniq > snp.ls
# 4048 SNPs with potential batch effects to be examed.

# on the nc2 server

nohup SNP_cluster_plot.R b01 plates1-53/     &
nohup SNP_cluster_plot.R b02 plates54-105/   &
nohup SNP_cluster_plot.R b03 plates106-156/  &
nohup SNP_cluster_plot.R b04 plates157-209/  &
nohup SNP_cluster_plot.R b05 plates210-261/  &
nohup SNP_cluster_plot.R b06 plates262-318/  &
nohup SNP_cluster_plot.R b07 plates319-367/  &

# It takes about 8~9 hours to classify the SNP callings of each batch,
# about 1 hour to grab sub-tables for 4000 SNPs,
# and 1.5 hours to plot 4000 SNPs.

# counting classes:
tail -n +2  b01/Ps.performance.txt | awk '{print $16}'  | sort | uniq -c
# and b02, b03 etc ...

#                           b01     b02     b03     b04     b05     b06     b07
# PolyHighResolution     521250  519578  518984  513339  516908  516373  515777
# NoMinorHom             104752  107367  109141  110955  102307  101447  102348
# MonoHighResolution      70838   66297   65634   67220   69333   76164   74700
# Hemizygous               1227    1227    1227    1227    1227    1227    1227
# OTV                      3437    4393    4279    3947    4218    4004    3901
# Other                   77710   80636   79786   82277   85380   79715   81393
# CallRateBelowThreshold   2723    2439    2886    2972    2564    3007    2591

# PolyHighResolution     66.2 +- 0.34 %
# NoMinorHom             13.5 +- 0.48 %
# MonoHighResolution      9.0 +- 0.53 %
# Hemizygous              0.2 +- 0    %
# OTV                     0.5 +- 0.04 %
# Other                  10.4 +- 0.31 %
# CallRateBelowThreshold  0.4 +- 0.03 %

# Use
SNP_class_pie.R
# to get a pie plot of different SNP classes.

#######################################
# For each SNP, in addition to its clustering plots, we also want a square plot
# to show us the calling classes across the 7 batches.
SNP_class_squares.R

#######################################
# Now, combine the square plots and the clustering plots.

IFS=$'\n'  snps=($(cat snp.ls))
# or in bash, it can also be
# readarray snps < ./snp.ls

for snp in ${snps[@]} ; do
    echo $snp
    convert ${snp}_class.png b01/$snp.png b02/$snp.png  +append r1.png
    convert b03/$snp.png b04/$snp.png b05/$snp.png      +append r2.png
    convert b06/$snp.png b07/$snp.png ${snp}_class.png  +append r3.png
    convert r1.png r2.png r3.png -append ${snp}_comb.png
done

mkdir batch_eff_png
mv *comb.png batch_eff_png

#######################################
# for the plate effect SNPs, do the same
cp snp.ls snp_ls.bak

tail -n +2 /kuser/shared/data/GWAS_backup/full_data/plate-effect/\
variant_plate_effects_v2.txt | awk '{print $3}' | sort | uniq > t.ls
# 30570 SNPs

grab -v -f snp_ls.bak t.ls > snp.ls
# 29903 SNPs not plotted yet

nohup SNP_cluster_plot.R b01 plates1-53/     &
nohup SNP_cluster_plot.R b02 plates54-105/   &
nohup SNP_cluster_plot.R b03 plates106-156/  &
nohup SNP_cluster_plot.R b04 plates157-209/  &
nohup SNP_cluster_plot.R b05 plates210-261/  &
nohup SNP_cluster_plot.R b06 plates262-318/  &
nohup SNP_cluster_plot.R b07 plates319-367/  &

SNP_class_squares.R


IFS=$'\n'  snps=($(cat snp.ls))

for snp in ${snps[@]} ; do
    echo $snp
    convert ${snp}_class.png b01/$snp.png b02/$snp.png  +append r1.png
    convert b03/$snp.png b04/$snp.png b05/$snp.png      +append r2.png
    convert b06/$snp.png b07/$snp.png ${snp}_class.png  +append r3.png
    convert r1.png r2.png r3.png -append ${snp}_comb.png
done

mkdir plate_eff_png
mv *_comb.png plate_eff_png
# about 12G of files
################################################################################
