#!/usr/bin/env python3

# tail -n +2 plink.genome | awk '{print $2,$4}' > t.in

edges = []
ifile =open('t.in')
for line in ifile:
    cols = line[:-1].split()
    edges.append([cols[0],cols[1]])
ifile.close()

import collections
node_dgres = collections.Counter()
nodes_1 = [x[0] for x in edges]
nodes_2 = [x[1] for x in edges]

node_dgres.update(nodes_1)
node_dgres.update(nodes_2)

# lets remove nodes according to their connection degrees
to_remove_list = []
for l in range(10000):
    if edges == []:
        break
    # find the most connected node
    to_remove_id = node_dgres.most_common(1)[0][0]
    to_remove_list.append(to_remove_id)

    # update edge set
    new_edges = [x for x in edges if to_remove_id not in x]
    edges = new_edges

    # update node connection degree
    node_dgres = collections.Counter()
    nodes_1 = [x[0] for x in edges]
    nodes_2 = [x[1] for x in edges]
    node_dgres.update(nodes_1)
    node_dgres.update(nodes_2)

for id in to_remove_list:
    print(id)