################################################################################
#                      CKB SNP Clustering (manual) check                       #
################################################################################

################################################################################
# Given batches of SNP calling from AxiomGT1, regressions have be performed to #
# detect SNPs with potential batch and other effects.                          #
# Some manual (visual) checking might be necessary for deciding thresholds     #
# and removal of SNPs with bad calling.                                        #
################################################################################

# on the nc2 server:
cd /kuser/kuangl/dev/ckb_batch_check

################################################################################
# 1. Input files                                                               #
################################################################################

################################################################################
# 1.1 SNP lists:

# plate effect:
/kuser/shared/data/GWAS_backup/full_data/plate-effect/\
variant_plate_effects_v2.txt
# batch   plate   probeset p-value
# 33621 entries, 30570 uniq SNPs

# batch effect:
/kuser/shared/data/GWAS_backup/full_data/batch_test/\
variant_batch_effects.txt
# batch   probeset        P-val
# 6407 entries, 4048 uniq SNPs
/kuser/shared/data/GWAS_backup/full_data/batch_test/\
variant_batch_norel_effects.txt
# batch   probeset        P-val
# 4154 entries, 2876 uniq SNPs

awk '{print $3}' /kuser/shared/data/GWAS_backup/full_data/plate-effect/\
variant_plate_effects_v2.txt    >  t.ls
awk '{print $2}' /kuser/shared/data/GWAS_backup/full_data/batch_test/\
variant_batch_effects.txt       >> t.ls
awk '{print $2}' /kuser/shared/data/GWAS_backup/full_data/batch_test/\
variant_batch_norel_effects.txt >> t.ls
sort t.ls | uniq | grep -v probeset > snp.ls
# 44183 entries, 34394 unique SNPs

################################################################################
# 1.2 calling files:

# on the nc2 server again:
# 7 batches of calling files at
/kuser/shared/data/GWAS_backup/

plates1-53/
plates54-105/
plates106-156/
plates157-209/
plates210-261/
plates262-318/
plates319-367/

# From each batch, we need the four files:
AxiomGT1.calls.txt
AxiomGT1.confidences.txt
AxiomGT1.snp-posteriors.txt
AxiomGT1.summary.txt

# and the plate cel file lists for plate highlighting
plates1-53.txt
plates54-105.txt
plates106-156.txt
plates157-209.txt
plates210-261.txt
plates262-318.txt
plates319-367.txt

################################################################################
# 2. SNP cluster plots                                                         #
################################################################################

################################################################################
# 2.1 SNPolisher SNP classification

mkdir b01 b02 b03 b04 b05 b06 b07

cat /kuser/shared/data/GWAS_backup/full_data/*stage1.bim | awk '{print $2'} | \
    sort | uniq > full_snp.ls
# 687236 SNPs

# install R package SNPolisher

SNP_classify.R b01 plates1-53/
# 1 calculate SNP clustering metrics for all SNPs in 'full_snp.ls'.
# 2 classify them into 7 categories.
# 3 grab the calls, confs, posterior and summary sub-tables for the listed SNPs
#    using SNPolisher's perl script.
# 4 plot (if uncommented the last two sections.)

# on the nc2 server
nohup SNP_classify.R b01 plates1-53/     &
nohup SNP_classify.R b02 plates54-105/   &
nohup SNP_classify.R b03 plates106-156/  &
nohup SNP_classify.R b04 plates157-209/  &
nohup SNP_classify.R b05 plates210-261/  &
nohup SNP_classify.R b06 plates262-318/  &
nohup SNP_classify.R b07 plates319-367/  &

# To get the 'metrics.txt' and 'Ps.performance.txt' files:
# It takes about 8~9 hours to classify all SNP callings of each batch.

# Then about 1 hour to grab sub-tables for 4000 SNPs,
# and 1.5 hours to plot 4000 SNPs.

# The plotting part is slow, and the pictures are hard to look at.

################################################################################
# 2.2 counting SNP calling classes

tail -n +2  b01/Ps.performance.txt | awk '{print $16}'  | sort | uniq -c
# and b02, b03 etc ...
for batch in b01 b02 b03 b04 b05 b06 b07 ; do
    tail -n +2 $batch/Ps.performance.txt | awk '{print $16}'  | sort | uniq -c
done

#                           b01     b02     b03     b04     b05     b06     b07
# PolyHighResolution     521250  519578  518984  513339  516908  516373  515777
# NoMinorHom             104752  107367  109141  110955  102307  101447  102348
# MonoHighResolution      70838   66297   65634   67220   69333   76164   74700
# Hemizygous               1227    1227    1227    1227    1227    1227    1227
# OTV                      3437    4393    4279    3947    4218    4004    3901
# Other                   77710   80636   79786   82277   85380   79715   81393
# CallRateBelowThreshold   2723    2439    2886    2972    2564    3007    2591

# PolyHighResolution      66.2 +- 0.34 %
# NoMinorHom              13.5 +- 0.48 %
# MonoHighResolution       9.0 +- 0.53 %
# Hemizygous               0.2 +- 0    %
# OTV                      0.5 +- 0.04 %
# Other                   10.4 +- 0.31 %
# CallRateBelowThreshold   0.4 +- 0.03 %

# Use
SNP_class_pie.R
# to get a pie plot of different SNP classes.

################################################################################
# 2.3 calling class square plots

# For each SNP, in addition to its clustering plots, we also want a square plot
# to show the calling classes across the 7 batches.

mkdir class_png

# to get SNPs classifications across all 7 batches into 't.in'
awk '{print $1,$16}' b01/Ps.performance.txt > t.out
grab -f snp.ls t.out > t.in

for batch in  b02 b03 b04 b05 b06 b07 ; do
    awk '{print $1,$16}' $batch/Ps.performance.txt > t.out
    grab -f snp.ls t.out > t2.in
    paste t.in t2.in > t.out

    # check SNP ids match
    awk '{if ($1 != $(NF-1) ) print "WAT?"}' t.out

    awk '{print $2}' t2.in > t3.in
    paste t.in t3.in > t2.in
    mv t2.in t.in
done

nohup SNP_class_squares.R &
# more than 200 SNPs per min
# 300,000 per day
# fast enough

################################################################################
# 2.4 re-visit the plotting problem

# We want to plot SNP clustering, but found the script SNP_classify.R
# which uses the SNPolisher library, was too slow (about 1 png per min).

# Using the ggplot2 library allows more freedom on manipulating the pictures.
# And it should be faster.

###################
# double check the orders of chips are the same between the
# call and summary files
for batch in b01 b02 b03 b04 b05 b06 b07 ; do
    echo $batch
    head -n 1 $batch/calls.txt > t.ls
    transpos_file t.ls > t1.ls
    head -n 1 $batch/summary.txt > t.ls
    transpos_file t.ls > t2.ls
    diff t1.ls t2.ls | wc
done
# no difference found, all are consistent.

###################
# double check the summary files has the order of A and B of the same SNP
# say, 'AX-100002645-A' followed by 'AX-100002645-B'
# then 'AX-100002667-A' followed by 'AX-100002667-B'

for batch in b01 b02 b03 b04 b05 b06 b07 ; do
    echo $batch
    awk '{print $1}' $batch/summary.txt | tail -n +2 > t.ls
    awk -F"-" '{print $2}' t.ls | uniq -c | awk '{print $1}' | uniq -c
    # all 2
    awk -F"-" '{print $3}' t.ls | sort | uniq -c
    # 'A' or 'B', nothing else
    awk -F"-" '{print $3}' t.ls | uniq -c | awk '{print $1}' | uniq -c
    # all 1, 'A' and 'B' are never consecutive
done

#######################################
# 2.4.1 get the sub-tables from the posterior files

for batch in b01 b02 b03 b04 b05 b06 b07 ; do
    get_posterior.py $batch > ${batch}.posterior
done

#######################################
# 2.4.2 generate avm files
# The avm (A_vs_M) files should have the a and b signals,
# in the form of A <- (log2(a) + log2(b))/2 and M <- log2(a) - log2(b).
# The calls, included into an avm file, are '0', '1', '2' and '-1' for missing.
# Here we change '-1' to '3' for missing.

# This script use the calls.txt and summary.txt files in the batch directories.
nohup get_avm.py b01 &
nohup get_avm.py b02 &
nohup get_avm.py b03 &
nohup get_avm.py b04 &
nohup get_avm.py b05 &
nohup get_avm.py b06 &
nohup get_avm.py b07 &

# 500 SNP avm files per min per job
# That's far better than 1 or 2 SNPs per min using the SNPolisher script.
# Finished making all 687236 x 7 avm files in 17 hours.
# However, DO NOT DO IT.

# 4.8 million extra files make the system very very slow.

# given the A_vs_M files
nohup SNP_cluster_plot.R b01 &
nohup SNP_cluster_plot.R b02 &
nohup SNP_cluster_plot.R b03 &
nohup SNP_cluster_plot.R b04 &
nohup SNP_cluster_plot.R b05 &
nohup SNP_cluster_plot.R b06 &
nohup SNP_cluster_plot.R b07 &

# 100,000 png files per job over the weekend.
# One job produces 26 pictures per minute, about 1 per two second.

# It might be the huge number of avm files slowing down the system.
# The normal speed is 70 pictures per minute per job on smaller sets.

# So, about 60 times faster. Good. And the pictures are nicer.

################################################################################
# 2.5 combine the square plot and the clustering plots.

IFS=$'\n'  snps=($(cat snp.ls))

for snp in ${snps[@]} ; do
    echo $snp
    convert class_png/$snp.png b01/$snp.png b02/$snp.png  +append r1.png
    convert b03/$snp.png b04/$snp.png b05/$snp.png        +append r2.png
    convert b06/$snp.png b07/$snp.png class_png/$snp.png  +append r3.png
    convert r1.png r2.png r3.png -append ${snp}_comb.png
done

mkdir to_exam_png.bak/
mv *comb.png to_exam_png.bak/
# tgz of the directory 'to_exam_png.bak/' is available as
# '/kuser/kuangl/backup/snp_clust_plot_manual_exam_png.tgz'

################################################################################
# 3. manual check of the clustering plots                                      #
################################################################################

#######################################
# preparing

mkdir to_exam_png
mkdir to_exam

split -l 800 snp.ls
mv x?? to_exam/
# got the lists of SNPs to be examed.
# Well, don't do too many at one time.
# check 800 SNPs in a batch

# manual_chk_res.table is the results of manual checking
# In the second column, '0' means the SNP failed the check, '1' means pass.

#######################################
# check 800 SNPs in a batch

awk '{print $1}' manual_chk_res.table | sort | uniq > examed.ls

rm  to_exam_png/*.png

awk '{print $1}' to_exam/xaa | sort > t1.ls
grab -f examed.ls -v t1.ls > to_exam.ls
# 'grep -w ' if 'grab' is not available.

awk '{print "cp to_exam_png.bak/" $1 "_comb.png to_exam_png/"}' to_exam.ls | sh

geeqie to_exam_png/
# 'gwenview' or whatever picture viewer if 'geeqie' is not available.
# now check the png pictures in the to_exam_png directory
# delete pictures of badly called SNPs.

# after two passes, collect the checking results, add to the result table.
ls to_exam_png | sed 's/_comb.png//' > t.in
get_chk_res.py >> manual_chk_res.table

################################################################################
# 4. checking the missing rates and MAF againt failing manual QC               #
################################################################################

ped_dir=/kuser/shared/data/GWAS_backup/full_data

plink --bfile $ped_dir/1-53_stage1     --missing --freq  --out b01
plink --bfile $ped_dir/54-105_stage1   --missing --freq  --out b02
plink --bfile $ped_dir/106-156_stage1  --missing --freq  --out b03
plink --bfile $ped_dir/157-209_stage1  --missing --freq  --out b04
plink --bfile $ped_dir/210-261_stage1  --missing --freq  --out b05
plink --bfile $ped_dir/262-318_stage1  --missing --freq  --out b06
plink --bfile $ped_dir/319-367_stage1  --missing --freq  --out b07

# get the sub_tables
for batch in b01 b02 b03 b04 b05 b06 b07 ; do
    grab -f snp.ls -c 2 $batch.lmiss | awk '{print $2,$5}' > t.lmiss
    sort_table -f snp.ls t.lmiss > t_$batch.lmiss
done

for batch in b01 b02 b03 b04 b05 b06 b07 ; do
    grab -f snp.ls -c 2 $batch.frq | awk '{print $2,$5}' > t.frq
    sort_table -f snp.ls t.frq > t_$batch.frq
done

# into manual_chk_res.table
# add batch_min_p batch_norel_min_p plate_min_p
# add batch_num_detected batch_norel_num_detected plate_num_detected

# add 7 missing-call-rates
# add 7 maf

# add max_miss and mean_miss
# add min_maf and mean_maf

# save manual_chk_res.xlsx into text file t.in
tail -n +2 t.in | awk '{print $2,$5,$9}' | grep -v NA > t2.in
get_plate_p_miss_hm.py > t3.in

plot_man_qc.R
# produces p1.png p2.png ... p9.png

################################################################################
# 5. plate effect plotting                                                     #
################################################################################

# modify the get_avm.py
nohup get_highlight_avm.py b01 &
nohup get_highlight_avm.py b02 &
nohup get_highlight_avm.py b03 &
nohup get_highlight_avm.py b04 &
nohup get_highlight_avm.py b05 &
nohup get_highlight_avm.py b06 &
nohup get_highlight_avm.py b07 &

# 33620 avm files to be generated
# done in 45 minutes
ls b0?/*avm > t.ls
split -l 1682 t.ls

# simplify the SNP_cluster_plot.R script
# SNP_highlight_cluster_plot.R

for ifile in x?? ; do
    nohup SNP_highlight_cluster_plot.R $ifile  &
done

# much faster... not posterior eclipses to plot this time.
# about 100 png files per minute per job
# all done in 22 minutes

rm x??

mkdir -p to_exam_png.bak/
# after proper backup
rm to_exam_png.bak/*.png

mv b0?/*.png to_exam_png.bak/
rm b0?/*.avm

################################################################################
# 6. plate effect plots manual check                                           #
################################################################################

mkdir to_exam_png
mkdir to_exam

tail -n +2 manual_chk_res.table| \
    awk '{if ($2==0) print $1}' > failed_clustering_QC_snp.ls

ls to_exam_png.bak > t.ls
awk -F"_" '{print $1,$2}' t.ls > t.in
grab -f failed_clustering_QC_snp.ls -v t.in | awk '{print $1 "_" $2}' > t.ls
# in the list of 33621 png files, 20490 left

random_shuffle_lines.py t.ls > t2.ls

split -l 800 t2.ls
mv x?? to_exam

#######################################

awk '{print $1}' plate_man_qc.table | sort | uniq > examed.ls
rm to_exam_png/*.png

grab -f examed.ls -v  to_exam/xas  > to_exam.ls

awk '{print "cp to_exam_png.bak/" $1 " to_exam_png/"}' to_exam.ls | sh
geeqie to_exam_png/
# after two passes

ls to_exam_png  > t.in
get_chk_res.py >> plate_man_qc.table

# the first 4000 SNP-plates, with the SNPs past the manual clustering QC,
# 57 failed plate clustering QC.
# 80 failed in 6400
# 1.25 %
# I tried 800 SNP-plates, with SNPs failed manual clustering QC,
# 460 out of 800 failed plate clustering Qc.

################################################################################