from Bio import SearchIO
import glob
import pprint
import os
pp = pprint.PrettyPrinter(indent=4)

for file in glob.glob("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/foldlibs_blast_db/*.bls"):
    blast_file = open(file)
    try:
        for result in SearchIO.parse(blast_file, 'blast-text'):
            for hit in result:
                for hsp in hit:
                    if hsp.evalue == 0.0 or hsp.evalue <= 1e-6:
                        # pp.pprint(hsp.evalue)
                        pp.pprint(hit.id)
                        #hh_file = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_overlap/"+hit.id+".hhm"
                        hh_file = "/scratch0/NOT_BACKED_UP/dbuchan/HHSearch/"+hit.id+".a3m"

                        et_tdb = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/et_overlap/"+hit.id+".tdb"
                        et_eig = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/et_overlap/"+hit.id+".eig"
                        gt_tdb = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/foldlib_overlap/"+hit.id[1:5]+hit.id[5].upper()+"0.tdb"
                        try:
                            os.remove(hh_file)
                        except:
                            print("hh file missing")
                        try:
                            os.remove(et_tdb)
                        except:
                            print("et tdb file missing")
                        try:
                            os.remove(et_eig)
                        except:
                            print("et eig file missing")
                        try:
                            os.remove(gt_tdb)
                        except:
                            print("tdb file missing")
    except:
        print("bad formatted output")

    # if found != 1:
    #     not_homologues.append(file)
