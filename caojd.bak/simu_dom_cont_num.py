#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

import sys
while '/share/python/lib64/python' in sys.path :
    sys.path.remove('/share/python/lib64/python') 
# the cluster got this wrong dir
from numpy import *

dom_ids=open(sys.argv[1]).read().split()

ifile=open('index/cath_s35.condef')
for line in ifile:
    c_line=ifile.readline().strip()
    
    dom_id = line.split()[0][1:]
    if dom_id not in dom_ids: continue 
    
    dom_len = int (line.split()[1])
    real_cont_num = int (line.split()[2])
    real_cont_num = real_cont_num * 2
    
    ############################################################################
    # simulate dom contact numbers 

    simulate_cont_nums=[]
    for iter in range(10000):
        cont_num=0
        for i in range(dom_len):
            res_cont_num=0
            r_num=random.uniform(0,1)
            if   r_num < 272118  /1925243 : res_cont_num = 0; 
            elif r_num < 504400  /1925243 : res_cont_num = 1;
            elif r_num < 763254  /1925243 : res_cont_num = 2;
            elif r_num < 1037957 /1925243 : res_cont_num = 3;
            elif r_num < 1299416 /1925243 : res_cont_num = 4;
            elif r_num < 1514757 /1925243 : res_cont_num = 5;
            elif r_num < 1687271 /1925243 : res_cont_num = 6;
            elif r_num < 1812034 /1925243 : res_cont_num = 7;
            elif r_num < 1885401 /1925243 : res_cont_num = 8;
            elif r_num < 1914925 /1925243 : res_cont_num = 9;
            elif r_num < 1923243 /1925243 : res_cont_num =10;
            elif r_num < 1924986 /1925243 : res_cont_num =11;
            elif r_num < 1925220 /1925243 : res_cont_num =12;
            elif r_num < 1925242 /1925243 : res_cont_num =13;
            else  : res_cont_num =14;
            cont_num += res_cont_num;
        simulate_cont_nums.append(cont_num)
    
    simulate_cont_nums=array(simulate_cont_nums)
    
    mean_cont_num=mean(simulate_cont_nums)
    std_cont_num=std(simulate_cont_nums)
    
    z_score=(real_cont_num - mean_cont_num )/std_cont_num;
    print(dom_id,z_score)

    
ifile.close()
