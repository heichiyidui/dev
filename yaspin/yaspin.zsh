################################################################################
#                               YASPIN                                         #
#                                                                              #
#      Yet Another Secondary structure Predictor usIng Neual network           #
################################################################################

################################################################################
#                                                                              #
# The method was published in Bioinformatics. 21(2):152-9. (Epub 2004 Sep 17). #
#                                                                              #
# It uses a hidden neural network model for protein secondary structure        #
# prediction.                                                                  #
#                                                                              #
# Three-states secondary structures are re-classified into 7 states:           #
# helix begins and ends, strand begins and ends are seperated from the helix   #
# and stand classes.                                                           #
#                                                                              #
# A simple back-propagation supervised feed-forward neural network was trained #
# for the prediction of emission probabilities of the 7 states using           #
# PSSMs (Position Specific Scoring Matrices) from PsiBlast searching.          #
#                                                                              #
# A hidden Markov model is then applied on the emission probabilities to       #
# determine the most likely transimissions between states.                     #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
# Neural network learning techniques improved a lot over the years.            #
# It will be interesting to check some new ideas on this old program.          #
#                                                                              #
# We have PSSMs, multiple alignments from Blast, DSSP and STRIDE definitions.  #
# We have CAO contact defintions.                                              #
# We want some model to predict secondary structure, residue exposure (contact #
# number) or even residue-residue side-chain contacts.                         #
#                                                                              #
################################################################################

################################################################################
# 1. The data set                                                              #
################################################################################

# CATH v4.1.0 S35
# 14421 domains
# 2204188 residues
# for QC steps, see ../cath/cath.zsh






