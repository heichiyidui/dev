################################################################################
# 1. apply the grm cutoff 0.025 on the old grm files                           #
################################################################################

################################################################################
# apply the cutoff 

for cohort in COR D1 D2 IR M3 NIT SLA UK; do
    gcta64 --grm grm_bak/${cohort}_0.5 \
           --grm-cutoff 0.025 \
           --make-grm --out grm/${cohort}
done

################################################################################
# REML estimation, with PCs

#PCs
#SLAGEN  4
#D1      0
#D2      4
#CORRIEL 0
#IRISH   1
#MGH     2
#NIT     1
#UK      1

########################################
# to get pc files 
####################
# 4 PCs 
for cohort in SLA D2 ; do 
    grep -v ID_1 sample/${cohort}_4PC.sample | grep -v "0 0 0 D C C C C B" | \
        awk '{print $2,$1,$5,$6,$7,$8}' > ${cohort}.qcovar
done
####################
# 2 PCs 
for cohort in M3; do 
    grep -v ID_1 sample/${cohort}_4PC.sample | grep -v "0 0 0 D C C C C B" | \
        awk '{print $2,$1,$5,$6}' > ${cohort}.qcovar
done
####################
# 1 PC 
for cohort in IR NIT UK; do 
    grep -v ID_1 sample/${cohort}_4PC.sample | grep -v "0 0 0 D C C C C B" | \
        awk '{print $2,$1,$5}' > ${cohort}.qcovar
done

mkdir pcs
mv *.qcovar pcs

########################################
# to get pheno files 
for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    grep -v ID_1 sample/${cohort}_4PC.sample | grep -v "0 0 0 D C C C C B" | \
        awk '{print $2,$1,$9}' > ${cohort}.phen
done

mkdir phen 
mv *.phen phen

########################################
# to get hsq 

prev=0.00003

for cohort in D2 IR M3 NIT SLA UK; do 
    gcta64 --grm grm/${cohort} \
        --pheno phen/${cohort}.phen --qcovar pcs/${cohort}.qcovar\
        --prevalence ${prev} \
        --reml --out ${cohort}
done

for cohort in D1 COR; do 
    gcta64 --grm grm/${cohort} --pheno phen/${cohort}.phen \
        --prevalence ${prev} \
        --reml --out ${cohort}
done


for cohort in COR D1 D2 IR M3 NIT SLA UK; do
    cp ${cohort}.hsq hsq/${cohort}_00003.hsq 
done 


################################################################################
# 2. get grm files from orginal bed files (no metatable filter etc.)           #
################################################################################

peddir=~ifoghbrc/scratch/CGTA-POLY/PLINK-FORMAT 

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    for chr in 01 02 03 04 05 06 07 08 09 10 \
               11 12 13 14 15 16 17 18 19 20 21 22; do 
        plink --noweb --file ${peddir}/${cohort}_IMP_chr_${chr} \
              --recode --transpose --out ${cohort}_chr_${chr}
    done
done

# got problems with COR cohort, missing is N instead of 0 

for cohort in COR; do 
    for chr in 01 02 03 04 05 06 07 08 09 10 \
               11 12 13 14 15 16 17 18 19 20 21 22; do 
        plink --noweb --file ${peddir}/${cohort}_IMP_chr_${chr} \
              --missing-genotype N \
              --recode --transpose --out ${cohort}_chr_${chr}
        sed -i 's/N N/0 0/g' ${cohort}_chr_${chr}.tped 
    done 
done 

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    cat ${cohort}_chr_??.tped > ${cohort}.tped &
    cp ${cohort}_chr_22.tfam ${cohort}.tfam
done

#######################################
# to remove duplications in SNP ids and bp positions from the tped files

# rmSNP_POS_dup.py:

#--------------------------------------
#!/usr/bin/env python2.7
import sys
 
ifile=open(sys.argv[1])
ids=set([])
bps=set([])
for line in ifile:
    cols=line.split()
    (id,bp)=(cols[1],cols[0]+'_'+cols[3]) # bp is actually chr_bp here
    if id in ids:
        continue;
    if bp in bps:
        continue 
    print(line[:-1])
    ids.add(id)
    bps.add(bp)
ifile.close()
#--------------------------------------

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    ./rmSNP_POS_dup.py ${cohort}.tped > temp_${cohort}.tped & 
done

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    mv temp_${cohort}.tped ${cohort}.tped
done

#######################################
# to get the plink binary files

mkdir bed_imp 

for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    plink --noweb --tfile ${cohort} --make-bed --out bed_imp/${cohort} &
done
for cohort in COR D1 D2 IR M3 NIT SLA UK; do 
    qsub t.sh $cohort 
done

