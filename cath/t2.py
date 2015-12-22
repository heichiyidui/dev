#!/usr/bin/env python3 
import sys 

import fastaio

seqs=fastaio.read_fasta(sys.argv[1])
print(seqs[0].split()[1])