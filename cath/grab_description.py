#!/usr/bin/env python3 

ifile=open('CathDomainDescriptionFile')

for line in ifile:
    if line.startswith('DOMAIN'):
        id=line.split()[1]
        name_str=''
        source_str=''
        for line in ifile:
            if line.startswith('NAME'):
                name_str+=line[10:-1]
            if line.startswith('SOURCE'):
                source_str+=line[10:-1]
            if line.startswith('//'):
                print(id,name_str,source_str)
                break;
ifile.close()
