#!/usr/bin/env python3

AA_CODE  = "ARNDCQEGHILKMFPSTWYV"
AA_2_int = {}
for i in range(len(AA_CODE)):
    AA_2_int[AA_CODE[i]] = i

CONT_GROUPS = \
    [['WA','WR','WN','WD','WQ','WE','WH','WI','WL','WK','WM','WS','WT','WV'],\
    ['FA','FR','FN','FD','FQ','FE','FG','FH','FK','FS','FT','WG','YA','YR',  \
     'YN','YD','YQ','YE','YG','YH','YK','YS','YT'],                          \
    ['II','IL','IM','IV','LI','LL','LM','LV','MI','ML','MM','MV','FI','FL',  \
     'FM','FV','YI','YL','YM','YV','VI','VL','VM','VV'],                     \
    ['AF','AW','AY','CW','IF','IW','IY','LF','LW','LY','MF','MW','MY','FF',  \
     'FW','FY','SF','SW','SY','TF','TW','TY','WF','WW','WY','YF','YW','YY',  \
     'VF','VW','VY'],                                                        \
    ['PA','PR','PN','PD','PC','PQ','PE','PG','PH','PI','PL','PK','PM','PF',  \
     'PS','PT','PW','PY','PV'],                                              \
    ['AP','RP','NP','DP','CP','QP','EP','GP','HP','IP','LP','KP','MP','FP',  \
     'PP','SP','TP','WP','YP','VP'],                                         \
    ['AQ','AE','IR','IQ','IE','IK','LR','LQ','LE','LK','MR','MQ','ME','MK',  \
     'SQ','SE','TQ','TE','VR','VQ','VE','VK'],                               \
    ['IA','IN','ID','IG','IS','IT','LA','LN','LD','LG','LS','LT','MA','MN',  \
     'MD','MG','MS','MT','VA','VN','VD','VG','VS','VT'],                     \
    ['AH','RH','NH','DH','CH','QH','EH','GH','IH','LH','KH','MH','SH','TH',  \
     'VH'],                                                                  \
    ['AR','AK','RA','RR','RN','RD','RQ','RE','RK','RS','RT','DR','DK','QA',  \
     'QR','QN','QQ','QK','QS','QT','ER','EQ','EK','KA','KR','KN','KD','KQ',  \
     'KE','KK','KS','KT','SR','SK','TR','TK'],                               \
    ['AN','AD','AG','RG','NN','ND','NG','DN','DD','DG','QD','QE','QG','EN',  \
     'ED','EE','EG','GN','KG','SN','SD','SG','TN','TD','TG'],                \
    ['HA','HR','HN','HD','HQ','HE','HG','HH','HI','HL','HK','HM','HF','HS',  \
     'HT','HW','HY','HV'],                                                   \
    ['RI','RL','RM','RV','QI','QL','QM','QV','KI','KL','KM','KV'],           \
    ['RF','RW','RY','DF','DW','DY','QF','QW','QY','EF','EW','EY','KF','KW',  \
     'KY'],                                                                  \
    ['NA','NR','NQ','NE','NI','NL','NK','NM','NF','NS','NT','NW','NY','NV'], \
    ['DA','DQ','DE','DI','DL','DM','DS','DT','DV','EA','EC','EI','EL','EM',  \
     'ES','ET','EV'],                                                        \
    ['AC','RC','NC','DC','CC','QC','GC','HC','IC','LC','KC','MC','FC','SC',  \
     'TC','WC','YC','VC'],                                                   \
    ['CA','CR','CN','CD','CQ','CE','CG','CI','CL','CK','CM','CF','CS','CT',  \
     'CY','CV'],                                                             \
    ['GA','GR','GD','GQ','GE','GG','GI','GL','GK','GM','GF','GS','GT','GW',  \
     'GY','GV'],                                                             \
    ['AA','AI','AL','AM','AS','AT','AV','SA','SI','SL','SM','SS','ST','SV',  \
     'TA','TI','TL','TM','TS','TT','TV']]

print('r1 r2 group')
for i in range(len(CONT_GROUPS)):
    for j in range(len(CONT_GROUPS[i])):
        print(CONT_GROUPS[i][j][0], CONT_GROUPS[i][j][1], i)

