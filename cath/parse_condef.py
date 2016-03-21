#!/usr/bin/env python3
import sys

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

for dom_id in domain_ids:
    #########################
    # read contact file
    contacts=[]
    ifile=open('conDef/'+dom_id)
    for line in ifile:
        contc=line.strip().replace(' ','-')
        contacts.append(contc)
    ifile.close()

    print(">"+dom_id)
    print(' '.join(contacts))
