#/usr/bin/env python

import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

ts=open(sys.argv[1], 'r').read()
#ts=open('t.ls', 'r').read()
seqlist=ts.split();

E_VALUE_THRESH = 0.001

for seqname in seqlist:
    fasta_string = open("seq/"+seqname+".fasta",'r').read()
    print fasta_string;
    result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string)
    blast_record = NCBIXML.read(result_handle)
    
    ofile=open("blastAln/"+seqname+".fasta",'w');
    ofile.write(fasta_string);
    
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                ofile.write(">1 \n"+hsp.query+'\n');
                ofile.write(">2 \n"+hsp.sbjct+'\n');
    
