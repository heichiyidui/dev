# scheme weights effect size estimated using inverse of the corresponding
# standard errors.

SCHEME STDERR

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER SNP
ALLELE A1 A2
EFFECT BETA
PVALUE P
STDERR SE
PROCESS st1.raw.assoc.linear
PROCESS st2.raw.assoc.linear
PROCESS st3.raw.assoc.linear
PROCESS st4.raw.assoc.linear
PROCESS st5.raw.assoc.linear


# === CARRY OUT AN INTERIM ANALYSIS OF THE FIRST FOUR FILES ===
OUTFILE ldl_direct .tbl
ANALYZE HETEROGENEITY

CLEAR

################################################################################
# SCHEME STDERR

# Load the metal results of direct LDL cohorts
MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT Effect
STDERR StdErr
PVALUE P-value

PROCESS ldl_direct1.tbl

# LOAD THE indirect LDL cohort
MARKER SNP
ALLELE A1 A2
EFFECT BETA
STDERR SE
PVALUE P

PROCESS st6.raw.assoc.linear

OUTFILE ldl_all .tbl
ANALYZE HETEROGENEITY

QUIT
