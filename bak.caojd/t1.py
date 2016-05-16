#!/usr/bin/env python3
#!/home/klinbrc/bin/python3

import sys

gap_cao=sys.argv[1]
end_cao=sys.argv[2]

import glob

ifile_names=glob.glob('x??_'+gap_cao+'_'+end_cao+'.out')

means=[]
for ifile_name in ifile_names:
    ifile=open(ifile_name)
    t_mean=float(ifile.read().strip())
    ifile.close()
    means.append(t_mean)

print(sum(means)/len(means))