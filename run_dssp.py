from Bio.PDB import *
import glob
from subprocess import call
import pprint
pp = pprint.PrettyPrinter(indent=4)

seqs = {}
# fasta= open("pdb_2015.fasta", "w")
for file in glob.glob("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_overlap/*.pdb"):
    pdb = file[-10:-4]

    call(["dssp", "-i", file, "-o", "tmp.dssp"])

    dssp_dict = make_dssp_dict("tmp.dssp")
    seq = ""
    ss_seq = ""
    print_ctl = 0
    for res in sorted(dssp_dict[0]):
        print(res)
        if 'A' in res[0]:
            print_ctl = 1
            res_info = (dssp_dict[0][res])
            seq = seq + res_info[0]
            ss_seq = ss_seq + res_info[1]
    #
    if print_ctl == 1 and len(seq) > 50:
        fasta_out = open("" + pdb + ".fasta", "w")
        fasta_out.write(">" + pdb + "\n")
        fasta_out.write(seq+"\n")
        fasta_out.close()

        ss_seq = ss_seq.replace("B", "C")
        ss_seq = ss_seq.replace("G", "C")
        ss_seq = ss_seq.replace("I", "C")
        ss_seq = ss_seq.replace("T", "C")
        ss_seq = ss_seq.replace("S", "C")
        ss_seq = ss_seq.replace("-", "C")

        dssp_out = open("" + pdb + "_A.ss", "w")
        dssp_out.write(seq+"\n")
        dssp_out.write(ss_seq)
        dssp_out.close()

    break
