import glob
from subprocess import call
import pprint
from multiprocessing import Pool
import shutil
pp = pprint.PrettyPrinter(indent=4)

def run_genth(file):
    pdb = file[-11:-6]
    #print(pdb)
    exe = "/scratch0/NOT_BACKED_UP/dbuchan/Applications/genthreader/GenThreader.sh"
    # # strsum_eigen 1jbeA.pdb 1jbeA.dssp $TDB_DIR/1jbeA.tdb $TDB_DIR/1jbeA.eig
    call([exe, "-i", file, "-j", pdb])
    for file in glob.glob(pdb+"*"):
        shutil.move(file, "dom_results")

seqs = {}
fasta_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/eigenthreader/seq_files/"
# fasta= open("pdb_2015.fasta", "w")
p = Pool(10)
p.map(run_genth, glob.glob(fasta_dir+"*.fasta"))
