import glob
import pprint
import csv
import re

pp = pprint.PrettyPrinter(indent=4)

omit_from_results_set = {}
omit_from_results_set_pdb = {}

with open("/home/dbuchan/eigendata/"
          "non_redundant_list_superfamily_level.txt") as non_redundant:
    for line in non_redundant:
        line = line.rstrip()
        omit_from_results_set[line] = 1
        omit_from_results_set_pdb[line[1:5]+line[5:6].upper()] = 1

# pp.pprint(omit_from_results_set_pdb)
# exit()

# get the list of scop folds
scop_list = {}
pdb_list = {}
with open("/home/dbuchan/eigendata/"
          "pdb_scop_class.txt") as scop_list_file:
    reader = csv.reader(scop_list_file, delimiter=',', quotechar='"')
    for row in reader:
        # print(row)
        scop_list[row[0]] = row[1]
        pdb = row[0][1:5]
        chain = row[0][5:6].upper()
        pdb_id = pdb+chain
        try:
            # pdb_list[pdb_id].append(r".".join(row[1].split(".")[:-2]))
            pdb_list.append(row[1])
        except:
            pdb_list[pdb_id] = [row[1], ]

# pp.pprint(scop_list)
# pp.pprint(pdb_list)
# exit()

bench_membership = {}
with open('/home/dbuchan/eigendata/'
          '/benchmark_family_members.txt') as benchfile:
    reader = csv.reader(benchfile, delimiter=',', quotechar='"')
    for row in reader:
        bench_membership[row[0]] = row[1]


def get_fifty_eigen(omit_from_results_set, scop_list, bench_membership):
    result_dir = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/results/"
    for file in glob.glob(result_dir+"optimised_t20_c9/*.out"):
        # if "1aba" not in file:
        #      continue
        print(file)
        results_list = []
        pdb = file[-9:-5]
        # print(pdb)
        chain = file[-5]
        pdb_id = pdb+chain
        # print(pdb_id)
        line_count = 0
        extra_hits = []
        if pdb_id in bench_membership:
            print(result_dir+"top_fifty_results/"+pdb_id+".top_fifty")
            out = open(result_dir+"top_fifty_results/"+pdb_id+".top_fifty",
                       "w+")
            out.write("# PDB ID: "+pdb_id+"\n")
            scop_class = bench_membership[pdb_id]
            out.write("# SCOP FAMILY: "+scop_class+"\n")
            # print(file)
            with open(file) as csvresult:
                lines = [line.split() for line in csvresult]
                lines.sort(key=lambda s: float(s[0]), reverse=True)
                # print(lines[len(lines)-5:len(lines)])
                for line in lines:
                    # print(line)
                    result_array = []
                    try:
                        scop_3_levels = ".".join(scop_class.split(".")[:-1])
                        this_3_levels = ".".join(scop_list[line[3]].split(".")[:-1])
                        # print(scop_3_levels)
                        scop_2_levels = ".".join(scop_class.split(".")[:-2])
                        this_2_levels = ".".join(scop_list[line[3]].split(".")[:-2])
                        # print(scop_2_levels)
                        # print(this_2_levels)
                        # print(fold_found)

                        if scop_list[line[3]] == scop_class:
                            pass
                            # print("FAMILY MATCH")
                        elif scop_3_levels == this_3_levels:
                            pass
                            # print("SUPERFAMILY MATCH")
                        elif scop_3_levels == this_3_levels:
                            print("FOUND HIT")
                            if line_count >= 50:
                                extra_hits.append(line)
                        else:
                            if line_count < 50:
                                out.write(line[0]+" "+line[1]+" "+line[2]+" "+line[3]+" "+scop_list[line[3]]+"\n")
                                line_count += 1

                    except:
                        print("Huh")
                for line in extra_hits:
                    out.write(line[0]+" "+line[1]+" "+line[2]+" "+line[3]+" "+scop_list[line[3]]+"\n")
        else:
            print("MISSING: "+pdb_id)
        #     # pp.pprint(results_list)
        #     results = results_list[len(results_list)-10:len(results_list)]
        #     # pp.pprint(results)
        #
        #     for element in reversed(results):
        #         out.write(element[0]+","+element[1]+","+element[2]+"\n")
        #         # pp.pprint(element)
        # else:
        #     continue
        # break



get_fifty_eigen(omit_from_results_set, scop_list, bench_membership)
