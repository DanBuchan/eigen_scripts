import re
import sys
import os.path
import glob
import shutil

overlapping_list = "/cs/research/bioinf/archive0/eigen_thread/joint_chains.txt"

with open(overlapping_list) as f:

    lines = f.read().splitlines()

chain_list = []
for line in lines:
    chain_list.append(line[0:4])

for file in glob.glob(r'/scratch0/NOT_BACKED_UP/dbuchan/HHSearch/*.a3m'):
    print(file)
    chain = file[42:46]
    print(chain)
    if chain in chain_list:
        shutil.copy(file, "/cs/research/bioinf/archive0/eigen_thread/HHSearch_a3m_overlap/")

#for file in glob.glob(r'/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/genthreader_chains/*.tdb'):
#    chain = file[81:85]
    #print(chain)
#    if chain in chain_list:
#        shutil.copy(file, "/cs/research/bioinf/archive0/eigen_thread/foldlib_overlap/")
