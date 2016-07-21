import glob
from subprocess import call
import pprint
from multiprocessing import Pool
import os

pp = pprint.PrettyPrinter(indent=4)

def run_strsum_eigen(file):
    pdb = file[-11:-4]
    scop_dir_id = file[-9:-7]
    # print(pdb)
    # print(scop_dir_id)
    scop_file = "/scratch0/NOT_BACKED_UP/dbuchan/pdbstyle-1.75/"+pdb+".ent"
    if os.path.isfile(scop_file):
        call(["dssp", "-i", scop_file, "-o", "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/dssp_data/"+pdb+".dssp"])
    else:
        print("MISSING: "+pdb)
    dssp = pdb+".dssp"
    tdb = pdb+".tdb"
    eig = pdb+".eig"
    dssp_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/dssp_data/"
    thread_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/et_data/"
    print(pdb+" > "+dssp)
    exe = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/eigenthreader/bin/strsum_eigen"
    # # strsum_eigen 1jbeA.pdb 1jbeA.dssp $TDB_DIR/1jbeA.tdb $TDB_DIR/1jbeA.eig
    call([exe, scop_file, dssp_dir+dssp, thread_dir+tdb, thread_dir+eig])

seqs = {}
pdb_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_hmm_complete/"
#run_strsum_eigen("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_overlap/d1fpoa2.hhm")
p = Pool(10)
p.map(run_strsum_eigen, glob.glob(pdb_dir+"*.hhm"))
