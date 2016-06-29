#!/usr/bin/env python3

ifile=open('t2.in')

p_range=[]
for p in range(0,100,5):
    p_range.append(p/1000000)

m_range=[]
for m in range(0,20,2):
    m_range.append(m/100)

boxes=[]
for i in range(20):
    boxes.append([0]*10)

pass_boxes=[]
for i in range(20):
    pass_boxes.append([0]*10)

ifile.readline()
for line in ifile:
    cols = line[:-1].split()
    pass_man_qc = int(cols[0])
    plate_p = float(cols[1])
    max_miss= float(cols[2])

    p_index = -1
    for p in p_range:
        if plate_p >= p: p_index +=1

    m_index = -1
    if max_miss >0.20: continue
    for m in m_range:
        if  max_miss >= m : m_index +=1

    boxes[p_index][m_index] += 1
    pass_boxes[p_index][m_index] += pass_man_qc
ifile.close()

print('plate_min_p max_miss pass_man_qc')
for p in range(20):
    for m in range(10):
        if boxes[p][m] ==0 :
            continue

        print(p*5/1000000,m*2/100,pass_boxes[p][m]/boxes[p][m])

