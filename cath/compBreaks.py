#!/usr/bin/env python3 
import sys

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])

ifile.close()

for did in domain_ids:
    ifile=open('cseq/'+did)
    cseq=''.join(ifile.read().split('\n')[1:])
    ifile.close()

    ifile=open('bseq/'+did)
    bseq=''.join(ifile.read().split('\n')[1:])
    ifile.close()

    
    aln_cseq = ''
    aln_bseq = ''
    j=0
    for i in range(len(cseq)):
        if cseq[i].isupper():
            if bseq[j].isupper():
                aln_cseq += cseq[i]
                aln_bseq += bseq[j]
                j += 1
            else: # bseq[j] == 'x':
                aln_cseq += '-'
                aln_cseq += cseq[i]
                aln_bseq += bseq[j]
                j += 1
                aln_bseq += bseq[j]
                j += 1
        else: # cseq[i].islower():
            if j >= len(bseq):
                aln_cseq += cseq[i]
                aln_bseq += '-'
            else:
                if bseq[j]=='x':
                    aln_cseq += cseq[i]
                    aln_bseq += bseq[j]
                    j += 1
                else:
                    aln_cseq += cseq[i]
                    aln_bseq += '-'
    
    for i in range(len(aln_bseq)):
        if aln_bseq[i] == 'x':
            print(aln_cseq[i])
    

    