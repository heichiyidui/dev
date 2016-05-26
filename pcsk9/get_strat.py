#!/usr/bin/env python3

s_ids = []
ifile=open('t.in')
for line in ifile:
    cols = line[:-1].split('\t')
    s_ids.append(cols[0])
ifile.close()
s_ids = set(s_ids)

#######################################
# let's read ages
age_base  = {}
age_resu1 = {}
age_resu2 = {}

for s_id in s_ids:
    age_base [s_id] = 'NA'
    age_resu1[s_id] = 'NA'
    age_resu2[s_id] = 'NA'

ifile=open('age_base.csv')
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    age_base[cols[0]] = cols[1]
ifile.close()

ifile=open('age_resu1.csv')
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    age_resu1[cols[0]] = cols[1]
ifile.close()

ifile=open('age_resu2.csv')
ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    age_resu2[cols[0]] = cols[1]
ifile.close()
# now we got all ages
#######################################
# read ldl-c

dir_ldl={}
ind_ldl={}
for s_id in s_ids:
    dir_ldl[s_id] = 'NA'
    ind_ldl[s_id] = 'NA'

ifile=open('direct_ldl_c.txt')
for line in ifile:
    cols = line[:-1].split()
    dir_ldl[cols[0]] = cols[1]
ifile.close()

ifile=open('indirect_ldl_c.txt')
for line in ifile:
    cols = line[:-1].split()
    ind_ldl[cols[0]] = cols[1]
ifile.close()

# now we got all ldl
#######################################
# decide strat

ifile=open('t.in')
for line in ifile:
    cols = line[:-1].split('\t')
    s_id      = cols[0]
    ascert    = cols[1]
    pass_qc   = cols[2]
    dir_base  = int(cols[3])
    dir_resu1 = int(cols[4])
    dir_resu2 = int(cols[5])
    in_dir    = int(cols[6])

    strat = 'NA'
    age   = 'NA'
    ldl   = 'NA'

    # failed QC
    if cols[2] == '0':
        cols.append(strat)
        cols.append(age)
        cols.append(ldl)
        print('\t'.join(cols))
        continue

    if dir_base + dir_resu1 + dir_resu2 > 0:
        strat = 'dir'
        ldl = dir_ldl [s_id]
        if   dir_base:
            age = age_base[s_id]
        elif dir_resu1:
            age = age_resu1[s_id]
        elif dir_resu2:
            age = age_resu2[s_id]
    elif in_dir:
        strat = '6'
        ldl = ind_ldl[s_id]
        age = age_resu2[s_id]

    if strat == 'dir':
        if   ascert == 'ICH':
            strat = '1'
        elif ascert == 'IS':
            strat = '2'
        elif ascert == 'SAH':
            strat = '3'
        elif ascert == 'MI/IHD':
            strat = '4'
        elif ascert == 'control':
            strat = '5'

    cols.append(strat)
    cols.append(age)
    cols.append(ldl)
    print('\t'.join(cols))
ifile.close()
