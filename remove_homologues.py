from Bio import SearchIO
import glob
import pprint
import os
pp = pprint.PrettyPrinter(indent=4)
hhcount = 0
etcount = 0
for file in glob.glob("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/foldlibs_blast_db/*.bls"):
    blast_file = open(file)
    try:
        for result in SearchIO.parse(blast_file, 'blast-text'):
            for hit in result:
                for hsp in hit:
                    if hsp.evalue == 0.0 or hsp.evalue <= 1e-6:
                        # pp.pprint(hsp.evalue)
                        pp.pprint(hit.id)
                        hh_file = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_hmm_reduced/"+hit.id+".hhm"
                        a3m_file = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_hmm_reduced/"+hit.id+".a3m"

                        et_tdb = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/et_reduced/"+hit.id+".tdb"
                        et_eig = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/et_reduced/"+hit.id+".eig"
                        try:
                            os.remove(hh_file)
                        except:
                            print("hh file missing")
                        try:
                            os.remove(a3m_file)
                        except:
                            print("a3m file missing")

                        try:
                            os.remove(et_tdb)
                        except:
                            print("et tdb file missing")
                        try:
                            os.remove(et_eig)
                        except:
                            print("et eig file missing")

    except:
        print("bad formatted output")

    # if found != 1:
    #     not_homologues.append(file)
print(hhcount)
print(etcount)
