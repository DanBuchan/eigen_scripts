import os
import glob
from itertools import product
from subprocess import call
import pprint

pp = pprint.PrettyPrinter(indent=4)

# export TDB_DIR=~/eigen_benchmark/et_overlap/
def run_eigen(file, vect):
    os.environ["TDB_DIR"] = "/home/dbuchan/eigen_benchmark/et_overlap/"
    pdb = file[-10:-4]
    print(pdb)
    seq_file = "/home/dbuchan/eigen_benchmark/seq_files/"+pdb+".fasta"
    ss2_file = "/home/dbuchan/eigen_benchmark/ss2_files/"+pdb+".ss2"
    out = "/home/dbuchan/eigen_benchmark/results/optimised/"+pdb+".out"
    stdout = "/home/dbuchan/eigen_benchmark/results/optimised/"+pdb+".stdout"
    stderr = "/home/dbuchan/eigen_benchmark/results/optimised/"+pdb+".stderr"

    eigen = "/home/dbuchan/eigen_benchmark/bin/eigenthreader"
    et_lst = "/home/dbuchan/eigen_benchmark/et.lst"
    f = open(stdout, "w")
    exe = "qsub"
    eigen_string = eigen+' -m -C0 -t20 -c9 -z1250 -F'+ss2_file+' '+seq_file+' '+file+' '+out+' '+et_lst
    #print(eigen_string)
    #call([exe, '-l', 'tmem=1.9G', '-l', 'h_vmem=1.9G', '-l', 'h_rt=24:0:0', '-o', stdout, '-e', stderr, '-b', 'y', '-v', 'TDB_DIR=/home/dbuchan/eigen_benchmark/et_overlap/', eigen, "-m"+pdb, "-C0", "-t20", "-c"+str(vect), "-z1250" "-F"+ss2_file, seq_file, file, out, et_lst])

seqs = {}
con_dir = "/home/dbuchan/eigen_benchmark/confiles/"
# fasta= open("pdb_2015.fasta", "w")

for vect in range(1, 21):
    for file in glob.glob(con_dir+"*.con"):
        run_eigen(file, vect)
        #break
    #break
        #p.starmap(run_eigen, zip( product(range(6, 21), glob.glob(con_dir+"*.con")) ) )
