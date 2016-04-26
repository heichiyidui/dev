#!/usr/bin/env python3
import sys
import math

BATCH_ID_MAP = {'1-53'   :'b01', '54-105' :'b02', '106-156':'b03', \
                '157-209':'b04', '210-261':'b05', '262-318':'b06', \
                '319-367':'b07'}

BATCH_ID = sys.argv[1]
# BATCH_ID = 'b01'
for b_id in BATCH_ID_MAP.keys():
    if BATCH_ID == BATCH_ID_MAP[b_id]:
        batch_long_id=b_id

# finding the set of cel files given a plate nor_id
nor_cel_map ={}
ifile=open('/kuser/shared/data/GWAS_backup/plates'+\
            batch_long_id+'.txt')
ifile.readline()
for line in ifile:
    nor_id = line[:-1].split('\\')[-2].split('_')[0]
    cel_id = line[:-1].split('\\')[-1]
    if nor_id not in nor_cel_map.keys():
        nor_cel_map[nor_id]=set([cel_id])
        continue
    nor_cel_map[nor_id].add(cel_id)
ifile.close()

# finding the plates to highlight given a SNP id
snp_to_check_plates={}
ifile=open('/kuser/shared/data/GWAS_backup/full_data/plate-effect/'+\
           'variant_plate_effects_v2.txt')
ifile.readline() # the header line
for line in ifile:
    cols = line[:-1].split()
    batch_long_id = cols[0][1:]
    batch_id = BATCH_ID_MAP[batch_long_id]
    if batch_id != BATCH_ID:
        continue
    nor_id   = cols[1]
    snp_id   = cols[2]
    if snp_id not in snp_to_check_plates.keys():
        snp_to_check_plates[snp_id]=[nor_id]
        continue
    snp_to_check_plates[snp_id].append(nor_id)
ifile.close()

################################################################################
# start generating avm files

in_call_file = open(BATCH_ID+'/calls.txt')
in_summ_file = open(BATCH_ID+'/summary.txt')

# read the header for cel_ids, highlight some of them later
cel_ids  = in_call_file.readline()[:-1].split()[1:]
in_summ_file.readline()

for line in in_call_file:
    cols = line[:-1].replace('\t-1','\t3').split()
    # want '3' for missing calls, not '-1'
    snp_id    = cols[0]
    snp_calls = [int(x) for x in cols[1:] ]

    line_a = in_summ_file.readline()[:-1]
    line_b = in_summ_file.readline()[:-1]

    if snp_id not in snp_to_check_plates.keys():
        continue

    log_a = [math.log(float(x),2) for x in line_a.split()[1:]]
    log_b = [math.log(float(x),2) for x in line_b.split()[1:]]

    for nor_id in  snp_to_check_plates[snp_id]:
        to_highlight_cel_ids = nor_cel_map[nor_id]

        ofile = open(BATCH_ID+'/'+snp_id+'_'+nor_id+'.avm','w')
        ofile.write('called\tM\tA\n')
        for i in range(len(snp_calls)):
            out_snp_call = str(snp_calls[i])
            if cel_ids[i] in to_highlight_cel_ids:
                out_snp_call = str(snp_calls[i]+4) # 4 5 6 7 instead of 0 1 2 3

            out_line = '\t'.join([out_snp_call, str(log_a[i]-log_b[i]), \
                                  str((log_a[i]+log_b[i])/2) ])
            ofile.write(out_line+'\n')
        ofile.close()

in_summ_file.close()
in_call_file.close()
################################################################################
