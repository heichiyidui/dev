################################################################################
# TREM2 analysis                                                               #
################################################################################

################################################################################
#                                                                              #
# First filter the probesets e.g. subset of probes that can be reliably        #
# detected in 80% of samples in at least one diagnostic group                  #
#                                                                              #
# Then, the effects of age, gender, collection site and RIN should be          #
# regressed out of the normalised gene expression data and the corresponding   #
# residuals compared between diagnostic groups (control, MCI and AD) using     #
# linear models followed by post-hoc t-tests.                                  #
#                                                                              #
################################################################################

########################################
# 1. subjects

# according to the 'Diagnosis at sampling TREM2+' column in 
# TREM2_cytokine_analysis_demographics.xlsx
# 64 AD
# 22 Control
# 17 MCI
#  1 FTD
#  3 mixed/DLB
#  3 PD
#  1 unknown

########################################
# 2. probes 

# 30 probes, 111 subjects, 3330 values
# 1392 undetected...

# need more than 51 AD detections or more than 17 CTL detections

