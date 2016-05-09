################################################################################
#                                                                              #
#                              C.A.O.J.D                                       #
#                                                                              #
################################################################################

################################################################################
# Scoring matrices for YASPIN and residue-residue contact prediction           #
################################################################################

################################################################################
# 1. CAO contact clustering

# To reduce the dimension of our problems, we need to further classify the
# residue contacts into groups?




################################################################################
# 2. select H class representative domains

mkdir index

# try get as many alignments as possible for each H class
wc ../cath/bl_out/* | sed 's/..\/cath\/bl_out\///' | grep -v total | \
    awk '{print $1/2 -1 ,$4}' | sort -r -g > t.in

sel_dom_repr.py > index/dom.ls
# 2099 domains
# 4428785 alignments
# 2110 alignments per domain, still smaller than mean of S35 domain hits.
# Small domain classes do have less Blast hits?

