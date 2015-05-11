#! /usr/bin/env python3
import struct

ifile=open('holl.ls')
ids=[]
for line in ifile:
    ids.append(line[:-1])
ifile.close()

isHoll=[]
ifile=open('../../PE/evokerSet/PE_affymetrix.sample')
ifile.readline()
ifile.readline()
numHoll=0;
for line in ifile:
    if line.split()[1] in ids:
        isHoll.append(True)
        isHoll.append(True)
        numHoll+=1
    else:
        isHoll.append(False)
        isHoll.append(False)
ifile.close();


for chr in ['01','02','03','04','05','06','07','08','09','0X','0Y','10','11',\
            '12','13','14','15','16','17','18','19','20','21','22','MT','XY']:

    ifile=open('../../PE/evokerSet/PE_'+chr+'_affymetrix.int.bin','rb')
    (rows,cols)=struct.unpack('ii', ifile.read(8))  
    
    ofile=open('PE_'+chr+'_affymetrix.int.bin','wb')
    ofile.write(struct.pack('ii',rows,numHoll*2))
    
    for i in range(rows):
        for j in range(cols):
            if isHoll[j]:
                ofile.write(ifile.read(4))
            else:
                ifile.read(4)
    ofile.close()
    ifile.close()

