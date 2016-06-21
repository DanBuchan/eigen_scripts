import pprint
import re
import sys
from collections import Counter

hh_list = "/cs/research/bioinf/archive0/eigen_thread/hhsearch_chain_list.txt"
gt_list = "/cs/research/bioinf/archive0/eigen_thread/psichain_sort.lst"

def get_list(this_list, path):
    with open(path) as infile:
        for line in infile:
            line = line.strip()
            if ".hhm" in line:
                if line[1:5] not in this_list:
                    this_list.append(line[1:5])
            elif ".tdb" in line:
                if line[0:4] not in this_list:
                    this_list.append(line[0:4])
    return(this_list)

hh_members = []
gt_members = []
hh_members = get_list(hh_members, hh_list)
gt_members = get_list(gt_members, gt_list)
chain_list = hh_members + gt_members
chains = Counter(chain_list)

for chain, count in chains.items():
    if count > 1:
        print(chain+","+str(count))
