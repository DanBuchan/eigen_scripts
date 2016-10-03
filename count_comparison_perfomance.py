import glob
import pprint
import csv
import re

pp = pprint.PrettyPrinter(indent=4)

omit_from_results_set = {}
omit_from_results_set_pdb = {}

with open("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/"
          "scop_data/non_redundant_list_superfamily_level.txt") as non_redundant:
    for line in non_redundant:
        line = line.rstrip()
        omit_from_results_set[line] = 1
        omit_from_results_set_pdb[line[1:5]+line[5:6].upper()] = 1

# pp.pprint(omit_from_results_set_pdb)
# exit()

# get the list of scop folds
scop_list = {}
pdb_list = {}
with open("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/"
          "scop_data/pdb_scop_class.txt") as scop_list_file:
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
with open('/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/'
          'scop_data/benchmark_family_members.txt') as benchfile:
    reader = csv.reader(benchfile, delimiter=',', quotechar='"')
    for row in reader:
        bench_membership[row[0]] = row[1]

# pp.pprint(bench_membership)
# exit()

# Parse eigenresults


def parse_eigen(omit_from_results_set, scop_list, bench_membership):
    result_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
                 "eigen_thread/results/"
    for file in glob.glob(result_dir+"optimised_t20_c9/*.out"):
        # if "1fvk" not in file:
        #     continue
        print(file)
        results_list = []
        pdb = file[-9:-5]
        print(pdb)
        chain = file[-5]
        pdb_id = pdb+chain
        print(pdb_id)
        if pdb_id in bench_membership:
            print(result_dir+"processed_comparison/"+pdb_id+".eigentop")
            out = open(result_dir+"processed_comparison/"+pdb_id+".eigentop",
                       "w+")
            out.write("# PDB ID: "+pdb_id+"\n")
            scop_class = bench_membership[pdb_id]
            out.write("# SCOP FAMILY: "+scop_class+"\n")
        #     # print(file)
            with open(file) as csvresult:
                lines = [line.split() for line in csvresult]
                lines.sort(key=lambda s: float(s[0]))
                # print(lines[len(lines)-5:len(lines)])
                for line in lines:
                    result_array = []
                    try:
                        scop_3_levels = ".".join(scop_class.split(".")[:-1])
                        this_3_levels = ".".join(scop_list[line[3]].split(".")[:-1])
                        # print(scop_3_levels)
                        if scop_list[line[3]] == scop_class:
                            print("FAMILY MATCH")
                        elif scop_3_levels == this_3_levels:
                            print("SUPERFAMILY MATCH")
                        else:
                            result_array = [line[0], line[3], scop_list[line[3]]]
                            results_list.append(result_array)
                    except:
                        result_array = [line[0], line[3], ""]
                        results_list.append(result_array)

            # pp.pprint(results_list)
            results = results_list[len(results_list)-10:len(results_list)]
            # pp.pprint(results)

            for element in reversed(results):
                out.write(element[0]+","+element[1]+","+element[2]+"\n")
                # pp.pprint(element)
        else:
            continue
        # break


def parse_hh(omit_from_results_set, scop_list, bench_membership):
    result_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0" \
                 "/eigen_thread/results/"
    start_parse_pattern = r"^\sNo\sHit"
    start_parse_re = re.compile(start_parse_pattern)
    end_parse_pattern = r"^No\s1"
    end_parse_re = re.compile(end_parse_pattern)
    length_pattern = r"^Match_columns\s(\d+)"
    length_re = re.compile(length_pattern)
    region_pattern = "(\d+)-(\d+)"
    region_re = re.compile(region_pattern)
    for file in glob.glob(result_dir+"hhresults/*.hhr"):
        seq_length = 0
        print(file)
        results_list = []
        pdb = file[-9:-5]
        print(pdb)
        chain = file[-5]
        pdb_id = pdb+chain
        print(pdb_id)
        parse_ctl = False
        if pdb_id in bench_membership:
            print(result_dir+"processed_comparison/"+pdb_id+".hhtop")
            out = open(result_dir+"processed_comparison/"+pdb_id+".hhtop", "w+")
            out.write("# PDB ID: "+pdb_id+"\n")
            scop_class = bench_membership[pdb_id]
            out.write("# SCOP FAMILY: "+bench_membership[pdb_id]+"\n")
            print(file)
            with open(file) as hhrResult:
                lines = hhrResult.read().splitlines()
                result_array = []
                for line in lines:
                    length_re_result = length_re.match(line)
                    start_parse_result = start_parse_re.match(line)
                    end_parse_result = end_parse_re.match(line)

                    if length_re_result:
                        seq_length = int(length_re_result.group(1))

                    if end_parse_result:
                        parse_ctl = False
                    if parse_ctl is True:
                        domain_id = line[4:12].lstrip().rstrip()
                        scop_family = line[12:21].lstrip()
                        prob = line[36:40].lstrip()
                        e_val = line[40:48].lstrip()
                        p_val = line[48:56].lstrip()
                        score = line[56:63].lstrip()
                        ss = line[63:70].lstrip()
                        cols = line[70:75].lstrip()
                        query_region = line[75:85].lstrip()
                        template_region = line[85:99].lstrip()

                        scop_3_levels = ".".join(scop_class.split(".")[:-1])
                        this_3_levels = ".".join(scop_family.split(".")[:-1])

                        overlap_size = 0
                        region_re_result = region_re.match(query_region)
                        if region_re_result:
                            overlap = int(region_result.group(2))-int(region_result.group(1))
                        # print(scop_3_levels)
                        if overlap_size/seq_length < 80:
                            continue

                        try:
                            # print(scop_3_levels)
                            if scop_family == scop_class:
                                print("FAMILY MATCH")
                            elif scop_3_levels == this_3_levels:
                                print("SUPERFAMILY MATCH")
                            else:
                                result_array = [score, domain_id, scop_family]
                                results_list.append(result_array)
                        except:
                            result_array = [score, domain_id, ""]
                            results_list.append(result_array)
                    if start_parse_result:
                        parse_ctl = True
                # pp.pprint(results_list)
                results = []
                if len(results_list) > 10:
                    results = results_list[0:10]
                else:
                    results = results_list
                # print(results)
                for element in results:
                    # pp.pprint(element)
                    element[2] = element[2].rstrip("(")
                    element[2] = element[2].rstrip()
                    out.write(element[0]+","+element[1]+","+element[2]+"\n")
                # pp.pprint(element)
        else:
            continue
        # break


def parse_genth(omit_from_results_set_pbd, pdb_list, bench_membership):
    result_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
                 "eigen_thread/results/"
    for file in glob.glob(result_dir+"genthreader_results/*.pgen.presults"):
        print(file)
        results_list = []
        pdb = file[-19:-15]
        print(pdb)
        chain = file[-15]
        pdb_id = pdb+chain
        print(pdb_id)
        if pdb_id in bench_membership:
            print(result_dir+"processed_comparison/"+pdb_id+".genthtop")
            out = open(result_dir+"processed_comparison/"+pdb_id+".genthtop",
                       "w+")
            out.write("# PDB ID: "+pdb_id+"\n")
            scop_class = bench_membership[pdb_id]
            # print(scop_class)
            out.write("# SCOP FAMILY: "+bench_membership[pdb_id]+"\n")
            # print(file)
            with open(file) as genthresult:
                lines = genthresult.read().splitlines()
                result_array = []
                for line in lines:
                    entries = line.split()
                    scop_3_levels = ".".join(scop_class.split(".")[:-1])
                    this_3_levels = ".".join(
                                        scop_list[entries[9]].split(".")[:-1])

                    # print(entries[9])
                    overlap_size = int(entries[7])-int(entries[6])
                    length = int(entries[8])
                    if overlap_size/length < 80:
                        continue

                    if scop_list[entries[9]] == scop_class:
                        print("FAMILY MATCH")
                    elif scop_3_levels == this_3_levels:
                        print("SUPERFAMILY MATCH")
                    else:
                        # print(entries[9][0:5])
                        try:
                            result_array = [entries[2], entries[9], scop_list[entries[9]]]
                            results_list.append(result_array)
                        except:
                            result_array = [entries[2], entries[9], ""]
                            results_list.append(result_array)

                results = results_list[0:10]

            for element in results:
                out.write(element[0]+","+element[1]+","+element[2]+"\n")
        #        pp.pprint(element)
        else:
            continue
        #break

parse_eigen(omit_from_results_set, scop_list, bench_membership)
parse_hh(omit_from_results_set, scop_list, bench_membership)
parse_genth(omit_from_results_set_pdb, pdb_list, bench_membership)
