#!/usr/bin/env python3 
import sys 
import random
random.seed(514)

cont_lines=open(sys.argv[1]).read().split('\n')
random.shuffle(cont_lines)

for line in cont_lines:
    if line != '':
        print(line)