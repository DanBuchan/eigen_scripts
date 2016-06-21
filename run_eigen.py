import glob
from subprocess import call
import pprint
from multiprocessing import Pool

pp = pprint.PrettyPrinter(indent=4)

def run_eigen(file):
    pdb = file[-10:-4]
    dssp = pdb+".dssp"
    print(pdb+" > "+dssp)
    call(["dssp", "-i", file, "-o", "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/dssp_overlap/"+dssp])


seqs = {}
pdb_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_overlap/"
# fasta= open("pdb_2015.fasta", "w")
p = Pool(10)
p.map_async(run_dssp, glob.glob(pdb_dir+"*.pdb"))
