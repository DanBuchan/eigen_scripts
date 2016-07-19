import glob
import pprint
import re
pp = pprint.PrettyPrinter(indent=4)

for file in glob.glob("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_hmm_complete/*.hhm"):
    pdb = file[-11:-4]
    #print(pdb)
    print_ctl = 0
    with open(file) as f:
        content = f.readlines()
        for line in content:
            line = line.rstrip()
            if line.startswith(">"):
                print_ctl = 0
            if line.startswith("#"):
                print_ctl = 0

            if line.startswith(">"+pdb):
                print_ctl = 1

            if print_ctl == 1:
                print(line)
    #break
