import csv
import pprint
import sys
import re
from Bio import SearchIO

pp = pprint.PrettyPrinter(indent=4)

scop_classification = {} # the SCOP class for each pdb domain ID
family_members = {} # A list of all pdbs in each SCOP class
superfamily_members = {} # A list of all pdbs in each SCOP fold
non_redundant_list = {} # all the benchmark members and their homologues at the
                        # SCOP family level (a.1.1.1)

non_redundant_list_superfamily = {} # all the benchmark members and their homologues at the
                                    # SCOP superfamily level (a.1)
bench_list = {}
# 'd9ximd_': 'g.3.1.1',
with open('/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/'
          'scop_data/dir.des.scop.1.75.txt') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
    for row in reader:
        if "px" in row[1]:
            domain_id = row[3]
            scop_classification[domain_id] = row[2]
            fold_id = ".".join(row[2].split(".")[:-1])
            if row[2] in family_members:
                family_members[row[2]].append(domain_id)
            else:
                family_members[row[2]] = [domain_id, ]

            if fold_id in superfamily_members:
                superfamily_members[fold_id].append(domain_id)
            else:
                superfamily_members[fold_id] = [domain_id, ]


csvfile.close()

# exit()

blast_results = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
                "eigen_thread/foldlibs_blast_db/"
b1 = open('./benchmark_family_members.txt', 'w+')
b2 = open('./benchmark_domain_id.txt', 'w+')
b3 = open('./benchmark_superfamily_members.txt', 'w+')

scop_pattern = "^(\w\.\d+\.\d+\.\d+)\s+"
scop_regex = re.compile(scop_pattern)
with open('/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/'
          'benchmark_names.txt') as names:
    for name in names:
        out_string = ''
        sup_out_string = ''
        bench_pdb = name[0:5]
        blast_file = blast_results+bench_pdb+".bls"
        bench_scop_family = ""
        try:
            for result in SearchIO.parse(blast_file, 'blast-text'):
                hit = result[0]
                hsp = hit[0]
                if hsp.evalue == 0.0 or hsp.evalue <= 1e-6:
                    result = scop_regex.match(hit.description)
                    bench_scop_family = result.group(1)
        except:
            print("bad formatted output")
        b2.write(bench_pdb+","+bench_scop_family+"\n")
        # b3.write(bench_pdb+","+".".join(bench_scop_family.split(".")[:-1])+"\n")

        non_redundant_list[bench_scop_family] = 1
        non_redundant_list_superfamily[bench_scop_family] = 1
        # print(bench_scop_family)
        if bench_scop_family in family_members:
            out_string += bench_pdb+","+bench_scop_family+","
            for dom_id in family_members[bench_scop_family]:
                non_redundant_list[dom_id] = 1
                out_string += dom_id+","
            out_string = out_string.rstrip(",")
            b1.write(out_string+"\n")
        else:
            print(bench_pdb+" NOT IN CLASSIFICATION OR NO HOMOLOGUES",
                  file=sys.stderr)

        bench_scop_superfamily = ".".join(bench_scop_family.split(".")[:-1])
        # print(bench_scop_superfamily)
        if bench_scop_superfamily in superfamily_members:
            sup_out_string += bench_pdb+","+bench_scop_superfamily+","
            for dom_id in superfamily_members[bench_scop_superfamily]:
                non_redundant_list_superfamily[dom_id] = 1
                sup_out_string += dom_id+","
            sup_out_string = sup_jout_string.rstrip(",")
            b3.write(sup_out_string+"\n")
        else:
            print(bench_pdb+" NOT IN CLASSIFICATION OR NO HOMOLOGUES",
                  file=sys.stderr)

        # break
names.close()

f1 = open('./pdb_scop_class.txt', 'w+')
for pdb in scop_classification:
    f1.write(pdb+","+scop_classification[pdb]+"\n")
f1.close()

f2 = open('./scop_class_members.txt', 'w+')
for family in family_members:
    f2.write(family)
    for pdb in family_members[family]:
        f2.write(","+pdb)
    f2.write("\n")
f2.close()

f6 = open('./scop_superfamily_members.txt', 'w+')
for family in superfamily_members:
    f6.write(family)
    for pdb in superfamily_members[family]:
        f6.write(","+pdb)
    f6.write("\n")
f6.close()

f3 = open("./non_redundant_list_family_level.txt", 'w+')
for pdb in non_redundant_list:
    f3.write(pdb+"\n")
f3.close()

f4 = open("./non_redundant_list_superfamily_level.txt", 'w+')
for pdb in non_redundant_list_superfamily:
    f4.write(pdb+"\n")
f4.close()
