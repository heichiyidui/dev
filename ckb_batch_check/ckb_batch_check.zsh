################################################################################
#                      CKB batch effect (manual) check                         #
################################################################################

################################################################################
# Given batches of SNP calling from AxiomGT1, regressions have be performed to #
# detect SNPs with potential batch and other effects.                          #
# Some manual (visual) checking might be neccessary for deciding thresholds    #
# and manual removal of SNPs with bad calling.                                 #
################################################################################

################################################################################
# 1. Input files

#######################################
# 1.1 SNP lists:

# on the nc2 server:
# plate effect:
/kuser/shared/data/GWAS_backup/full_data/\
plate-effect/variant_plate_effects_v2.txt
# 33621 entries, 30570 uniq SNPs

# batch effect:
/kuser/shared/data/GWAS_backup/full_data/\
batch_test/variant_batch_effects.txt
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

# with each batch, we need four files:
AxiomGT1.calls.txt
AxiomGT1.confidences.txt
AxiomGT1.snp-posteriors.txt
AxiomGT1.summary.txt

################################################################################

