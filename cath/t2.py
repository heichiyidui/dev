#!/usr/bin/env python3 
import sys 

ifile=open(sys.argv[1])

for line in ifile:
    if line.startswith('# Query:'):
        dom_id =  line.split()[2]
        ifile.readline()
        ifile.readline()
        ifile.readline()
        line = ifile.readline()
        no_hits = line.split()[1]
        print(dom_id,no_hits)
    
ifile.close()
