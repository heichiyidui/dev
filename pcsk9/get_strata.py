#!/usr/bin/env python3

################################################################################
# 1. read the ck_id list and sex from the fam file                             #
# 2. read the study id map
# 3. read the ldl files
# 4. read the age at test file
# 5. read the sample summary file
# 6. select, grouping and output

################################################################################
# 1. read FID and IID from the fam file.
ifile=open('ckb_ph12_s3.fam')
FID={}
IID=[]
SEX={}
for line in ifile:
    cols = line[:-1].split()
    IID.append(cols[1])
    FID[cols[1]]=cols[0]
    SEX[cols[1]]=cols[4]
ifile.close()
IID_set=set(IID)

################################################################################
# 2. read the study_id map and ascertainment map

#######################################
# Updating the ascertainment file, takes some time.
# Don't have to do it again after updating.
# ifile=open('GWAS_SNPdata_samples.csv')
# for line in ifile:
#     cols = line[:-1].split('\t')
#     iid_s = cols[0].split('-')[0] # base iid
#     is_iid_s_found =False
#     mapping_iid = 'NA' # the mapping IID in the fam file
#     for iid in IID_set:
#         if iid.startswith(iid_s):
#             is_iid_s_found = True
#             mapping_iid = iid
#             break;
#     if not is_iid_s_found:
#         continue
#     cols[0]=mapping_iid
#     print('\t'.join(cols))
# ifile.close()

iid_2_sid={}
sid_2_iid={}
ascert={}
for iid in IID:
    iid_2_sid[iid]='NA'
    ascert[iid]='NA'

ifile=open('GWAS_SNPdata_samples.csv')
for line in ifile:
    cols = line[:-1].split('\t')
    if cols[0] not in IID_set:
        continue
    iid_2_sid[cols[0]] = cols[1]
    sid_2_iid[cols[1]] = cols[0]
    ascert   [cols[0]] = cols[2]
ifile.close()

################################################################################
# 3. read the ldl files
dir_ldl={}
dir_ldl_survey={}
for iid in IID:
    dir_ldl[iid] = 'NA'
    dir_ldl_survey[iid]='NA' # 0, 1, 2 for baseline, resurvey 1 and resurvey 2

ifile=open('LDL-c_biochem_data.csv')
for line in ifile:
    cols = line[:-1].split('\t')
    is_iid_found = False
    if cols[0] in sid_2_iid.keys():
        iid = sid_2_iid[cols[0]]
        is_iid_found = True
    elif cols[1] in sid_2_iid.keys():
        iid = sid_2_iid[cols[1]]
        is_iid_found = True

    if not is_iid_found:
        continue
    if cols[6] == '0': # the 'use?' column
        continue
    if cols[3] == '': # the 'ldl_result_x10000' column
        continue
    dir_ldl_survey[iid] = cols[2]
    dir_ldl[iid] = cols[3]
ifile.close()

ind_ldl={}
ind_ldl_survey={}
for iid in IID:
    ind_ldl[iid]='NA'
    ind_ldl_survey[iid]='NA'
ifile=open('ind_ldl.csv')
for line in ifile:
    cols = line[:-1].split('\t')
    if cols[0] not in sid_2_iid.keys():
        continue
    iid = sid_2_iid[cols[0]]
    if cols[1] == '':
        continue
    ldl = str(float(cols[1])/100)
    ind_ldl[iid]=ldl
    ind_ldl_survey[iid]='2'
ifile.close()

################################################################################
# 4. read the age at test file
age_base  = {}
age_resu1 = {}
age_resu2 = {}
for iid in IID:
    age_base [iid] = 'NA'
    age_resu1[iid] = 'NA'
    age_resu2[iid] = 'NA'

ifile=open('age_base.csv')
for line in ifile:
    cols = line[:-1].split('\t')
    if cols[0] not in sid_2_iid.keys():
        continue
    iid = sid_2_iid[cols[0]]
    if cols[1] == '':
        continue
    age_base[iid]=cols[1]
ifile.close()

ifile=open('age_resu1.csv')
for line in ifile:
    cols = line[:-1].split('\t')
    if cols[0] not in sid_2_iid.keys():
        continue
    iid = sid_2_iid[cols[0]]
    if cols[1] == '':
        continue
    age_resu1[iid] = cols[1]
ifile.close()

ifile=open('age_resu2.csv')
for line in ifile:
    cols = line[:-1].split('\t')
    if cols[0] not in sid_2_iid.keys():
        continue
    iid = sid_2_iid[cols[0]]
    if cols[1] == '':
        continue
    age_resu2[iid] = cols[1]
ifile.close()
# Age of resurvey 2 for subject 880310350 is missing.
# Given the baseline age 50 and resurvey 1 age 52, put 58 into age_resu2.csv
# The most popular age gap combinations are 2,6 and 3,5.
################################################################################
# 5. read the sample summary file
# save the sheet 2 of PCSK9_sample_summary.xlsx into PCSK9_sample_summary.csv
# after deleting the column X

strata={}

strata['COPD case']     = 'NA'
strata['fatalIHD case'] = '4'
strata['ICH case']      = '1'
strata['ICH control']   = '5'
strata['IS case']       = '2'
strata['MI case']       = '4'
strata['Resurvey 2']    = '6'
strata['SAH case']      = '3'
strata['Unselected']    = '6'

ifile=open('PCSK9_sample_summary.csv')
ifile.readline()
ifile.readline()
print('FID IID sid rc stratum age sex ldl_c')

for line in ifile:
    cols = line[:-1].split('\t')
    sid = cols[0]
    if sid not in sid_2_iid.keys():
        continue
    iid = sid_2_iid[sid]
    # print('\t'.join([cols[1],ascert[sid_2_iid[sid]]]))
    # the first 'IS' should be 'ICH case' in PCSK9_sample_summary
    # otherwise perfect match to what we got already.

    if cols[2] == '1':
        continue # the excl column
    if cols[4] == '0':
        continue # the 'still OK' column

    stratum = 'NA'
    age='NA'
    ldl='NA'

    is_dir_found =False

    if '1' in cols[5:10]:
        is_dir_found = True
    if dir_ldl[iid] == 'NA':
        is_dir_found = False

    if is_dir_found :
        stratum = strata[ascert[iid]]
        ldl = dir_ldl[iid]
        if '0' == dir_ldl_survey[iid]:
            age = age_base[iid]
        elif '1' == dir_ldl_survey[iid]:
            age = age_resu1[iid]
        elif '2' == dir_ldl_survey[iid]:
            age = age_resu2[iid]

    if is_dir_found:
        print(FID[iid], iid, sid, sid[:2], stratum, age, SEX[iid], ldl)
        continue

    is_ind_found = False
    if '1' in cols[10:12]:
        is_ind_found = True
    if ind_ldl[iid] == 'NA':
        is_ind_found = False

    if is_ind_found:
        stratum = '6'
        ldl = ind_ldl[iid]
        age = age_resu2[iid]
        print(FID[iid], iid, sid, sid[:2], stratum, age, SEX[iid], ldl)
ifile.close()

################################################################################


