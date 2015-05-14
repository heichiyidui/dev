#!/usr/bin/env python3

from snn import Snn
net=Snn(3,2,2,0.01)
net.read_net('t.net')

import random
random.seed(514)

for i in range(1000):
    input=(random.normalvariate(0,1),\
           random.normalvariate(0,1),\
           random.normalvariate(0,1))
    target=net.propagate(input)
    print(input,target)
