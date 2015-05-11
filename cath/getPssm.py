#!/usr/bin/env python3
##!/home/klinbrc/bin/python3

import os
import sys
import subprocess

################################################################################
# global parameters                                                            #
################################################################################

#blast='~/bin/psiblast ' 

blast='/usr/bin/psiblast ' 

db=' -db ../nr/nrfilt '
dbtype=' -dbtype prot '

################################################################################
# 1. make the output directions                                                #
################################################################################

# the PSSM files from Blast searching
if not os.path.exists('./pssm'):
    subprocess.call('mkdir ./pssm',shell=True)

################################################################################
# 2. read the list of query sequence files                                     #
################################################################################

ifile=open(sys.argv[1]) 
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])
ifile.close()

################################################################################
# 3. PsiBlast the consensus sequences                                          #
################################################################################

for id in domain_ids:
    if not os.path.isfile('./cseq/'+id):
        sys.stderr.write('Failed to open the query sequence file seq/'+id+'\n')
        exit(1)
        
    print('PsiBlast searching...')
    subprocess.call(blast+db+' -num_iterations 3'+\
                    ' -query cseq/'+id+\
                    ' -out_ascii_pssm pssm/'+id,shell=True)
    
    
