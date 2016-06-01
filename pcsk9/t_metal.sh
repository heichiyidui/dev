# scheme weights effect size estimated using inverse of the corresponding
# standard errors.

SCHEME STDERR

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER SNP
ALLELE A1 A2
EFFECT BETA
PVALUE P
STDERR SE
PROCESS st1.assoc.linear
PROCESS st2.assoc.linear
PROCESS st3.assoc.linear
PROCESS st4.assoc.linear
PROCESS st5.assoc.linear


# === CARRY OUT AN INTERIM ANALYSIS OF THE FIRST FOUR FILES ===
OUTFILE pcsk9_direct .tbl
ANALYZE HETEROGENEITY


CLEAR

################################################################################
SCHEME STDERR

# Load the metal results of direct LDL cohorts
MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT Effect
STDERR StdErr
PVALUE P-value

PROCESS pcsk9_direct1.tbl

# LOAD THE indirect LDL cohort
MARKER SNP
ALLELE A1 A2
EFFECT BETA
STDERR SE
PVALUE P

PROCESS st6.assoc.linear



OUTFILE pcsk9_all .tbl
ANALYZE HETEROGENEITY

QUIT
