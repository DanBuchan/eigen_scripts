import os
import glob
from itertools import product
from subprocess import call
import pprint
from multiprocessing import Pool

pp = pprint.PrettyPrinter(indent=4)


def run_eigen(file):
    os.environ["TDB_DIR"] = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/et_overlap/"
    pdb = file[-11:-6]
    #print(pdb)
    seq_file = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/seq_files/"+pdb+".fasta"
    ss2_file = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/ss2_files/"+pdb+".ss2"
    con_file = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/confiles/"+pdb+".con"
    out = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/results/optimised_t20_c9/"+pdb+".out"
    stdout = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/results/optimised_t20_c9/"+pdb+".stdout"
    #
    exe = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/bin/eigenthreader"
    et_lst = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/et.lst"
    f = open(stdout, "w")
    print(exe+" -m"+pdb+" -c9"+" -C0"+" -t20"+" -z1250"+" -F"+ss2_file+" "+file+" "+con_file+" "+out+" "+et_lst)
    call([exe, "-m"+pdb, "-c9", "-C0", "-t20", "-z1250", "-F"+ss2_file, file, con_file, out, et_lst], stdout=f)

fasta_dir = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/seq_files/"
p = Pool(10)
p.map(run_eigen, glob.glob(fasta_dir+"*.fasta"))
