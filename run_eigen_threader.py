import os
import glob
from itertools import product
from subprocess import call
import pprint
from multiprocessing import Pool

pp = pprint.PrettyPrinter(indent=4)


def run_eigen(file):
    os.environ["TDB_DIR"] = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/et_reduced/"
    pdb = file[-11:-6]
    #print(pdb)
    seq_file = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/seq_files/"+pdb+".fasta"
    ss2_file = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/ss2_files/"+pdb+".ss2"
    con_file = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/contact_subsets/random/L/"+pdb+"_L.con"
    out = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/results/random/L/"+pdb+".out"
    stdout = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/results/random/L/"+pdb+".stdout"
    #
    exe = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/bin/eigenthreader"
    et_lst = "/scratch0/NOT_BACKED_UP/dbuchan/Applications/genthreader/data/et.lst"
    f = open(stdout, "w")
    print(exe+" -m"+pdb+" -c9"+" -C0"+" -t20"+" -z1250"+" -F"+ss2_file+" "+file+" "+con_file+" "+out+" "+et_lst)
    call([exe, "-m"+pdb, "-c9", "-C0", "-t20", "-z1250", "-F"+ss2_file, file, con_file, out, et_lst], stdout=f)

fasta_dir = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/seq_files/"
p = Pool(10)
p.map(run_eigen, glob.glob(fasta_dir+"*.fasta"))
