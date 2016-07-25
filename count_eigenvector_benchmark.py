import glob
import pprint
import csv
from operator import itemgetter

pp = pprint.PrettyPrinter(indent=4)

# load in the list of scop family members of the benchmarks set
omit_from_results_set = {}
with open("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/scop_data/non_redundant_list_family_level.txt") as non_redundant:
    for line in non_redundant:
        line = line.rstrip()
        omit_from_results_set[line] = 1

# get the list of scop folds
scop_list = {}
with open("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/scop_data/pdb_scop_class.txt") as scop_list_file:
    reader = csv.reader(scop_list_file, delimiter=',', quotechar='"')
    for row in reader:
        scop_list[row[0]] = row[1]

# pp.pprint(scop_list)

bench_membership = {}
with open('/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/scop_data/benchmark_family_members.txt') as benchfile:
    reader = csv.reader(benchfile, delimiter=',', quotechar='"')
    for row in reader:
        bench_membership[row[0]] = row[1]
# read in benchmark family members to find out which


for result_dir in glob.glob("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/results/distance/c*"):
    print(result_dir)
    for file in glob.glob(result_dir+"/*.out"):
        results_list = []
        pdb = file[-9:-5]
        chain = file[-5]
        pdb_id = pdb+chain
        #print(pdb_id)
        if pdb_id in bench_membership:
            # print(result_dir+"/"+pdb_id+".top")
            out = open(result_dir+"/"+pdb_id+".top", "w+")
            out.write("# PDB ID: "+pdb_id+"\n")
            out.write("# SCOP FAMILY: "+bench_membership[pdb_id]+"\n")
            # print(file)
            with open(file) as csvresult:
                lines = [line.split() for line in csvresult]
                lines.sort(key=lambda s: float(s[0]))
                # print(lines[len(lines)-5:len(lines)])

                for line in lines:
                    if line[3] in omit_from_results_set:
                        # print("SKIPPING")
                        pass
                    else:
                        try:
                            result_array = [line[0], line[3], scop_list[line[3]]]
                        except:
                            result_array = [line[0], line[3], ""]
                        results_list.append(result_array)
            results = results_list[len(results_list)-10:len(results_list)]
            for element in reversed(results):
                out.write(element[0]+","+element[1]+","+element[2]+"\n")
                # pp.pprint(element)
        else:
            continue
        # exit()
        # NOW process top 1 and top 2
        # pp.pprint(results_list[0:5])
