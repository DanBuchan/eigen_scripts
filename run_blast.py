import glob
from subprocess import call
import pprint
from multiprocessing import Pool

pp = pprint.PrettyPrinter(indent=4)

def run_blast(file):
    pdb = file[-11:-6]
    exe = "/scratch0/NOT_BACKED_UP/dbuchan/Applications/ncbi-blast-2.2.31+/bin/blastp"
    db = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/foldlibs_blast_db/hhsearch.fasta"
    # # strsum_eigen 1jbeA.pdb 1jbeA.dssp $TDB_DIR/1jbeA.tdb $TDB_DIR/1jbeA.eig
    call([exe, "-query", file, "-out", pdb+".bls", "-db", db])

seqs = {}
fasta_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/eigenthreader/seq_files/"
# fasta= open("pdb_2015.fasta", "w")
p = Pool(10)
p.map(run_blast, glob.glob(fasta_dir+"*.fasta"))
