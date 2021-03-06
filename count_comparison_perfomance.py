import glob
import pprint
import csv
import re

# family_removal_set = ['1i4jA', '1j3aA', '1chdA', '1dmgA', '1mk0A', '1wjxA',
#                       '1jkxA', '1lm4A', '1g2rA', '1k6kA', '1fcyA', '3dqgA',
#                       '1xkrA', '1dixA', '1nrvA', '1hxnA', '1dbxA', '1whiA',
#                       '1im5A', '1kq6A']
# superfamily_removal_set = ['1ej0A', '1ne2A', '1i4jA', '1gbsA', '5ptpA',
#                            '1xdzA', '1a3aA', '1k7jA', '1j3aA', '1fk5A',
#                            '1aoeA', '1behA', '1kqrA', '1i58A', '1chdA',
#                            '3borA', '1iwdA', '1dmgA', '1mk0A', '1npsA',
#                            '1gzcA', '1htwA', '1wjxA', '1ktgA', '1i1nA',
#                            '1jkxA', '1g9oA', '2hs1A', '1jbkA', '1ql0A',
#                            '1lm4A', '1g2rA', '1atzA', '1f6bA', '1gz2A',
#                            '1i71A', '1k6kA', '1bebA', '1hdoA', '1fcyA',
#                            '3dqgA', '1fqtA', '1dsxA', '1jyhA', '1h2eA',
#                            '1cxyA', '1xkrA', '1vhuA', '1d4oA', '1dixA',
#                            '1nrvA', '1vfyA', '1lo7A', '1hxnA', '1ckeA',
#                            '1cjwA', '1c44A', '1wkcA', '1qf9A', '1dbxA',
#                            '1tzvA', '1h0pA', '1c52A', '1aapA', '1ctfA',
#                            '1gmxA', '1bsgA', '1ek0A', '1ihzA', '1mugA',
#                            '1whiA', '1im5A', '1kq6A', '1tqhA']
pp = pprint.PrettyPrinter(indent=4)

def readOmitList(file):
    omit_from_results_set = {}
    omit_from_results_set_pdb = {}
    with open(file) as non_redundant:
        for line in non_redundant:
            line = line.rstrip()
            omit_from_results_set[line] = 1
            omit_from_results_set_pdb[line[1:5]+line[5:6].upper()] = 1
    return(omit_from_results_set, omit_from_results_set_pdb)


(omit_from_results_set,
 omit_from_results_set_pdb) = readOmitList("/mnt/bioinf/archive0/"
                                                  "eigen_thread/scop_data/"
                                                  "non_redundant_list_family"
                                                  "_level.txt")
(superfamily_omit_from_results_set,
 superfamily_omit_from_results_set_pdb) = readOmitList("/mnt/bioinf/archive0/"
                                                       "eigen_thread/"
                                                       "scop_data/"
                                                       "non_redundant_list_"
                                                       "superfamily_level.txt")
# pp.pprint(omit_from_results_set_pdb)
# exit()

# get the list of scop folds
scop_list = {}
pdb_list = {}
with open("/mnt/bioinf/archive0/eigen_thread/"
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
with open('/mnt/bioinf/archive0/eigen_thread/'
          'scop_data/benchmark_family_members.txt') as benchfile:
    reader = csv.reader(benchfile, delimiter=',', quotechar='"')
    for row in reader:
        bench_membership[row[0]] = row[1]

# pp.pprint(bench_membership)

def parse_eigen(omit_from_results_set, scop_list, bench_membership):
    result_dir = "/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/results/"
    for file in glob.glob(result_dir+"optimised_t20_c9/*.out"):
    #for file in glob.glob(result_dir+"c1/*.out"):
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
            print(result_dir+"processed_comparison_family/"+pdb_id+".eigentop")
            out = open(result_dir+"processed_comparison_family/"+pdb_id+".eigentop",
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
                    if float(line[4]) < 0.4:
                        continue

                    try:
                        scop_3_levels = ".".join(scop_class.split(".")[:-1])
                        this_3_levels = ".".join(scop_list[line[6]].split(".")[:-1])
                        # print(scop_3_levels)
                        if scop_list[line[6]] == scop_class:
                            print("FAMILY MATCH")
                        # elif scop_3_levels == this_3_levels:
                            # print("SUPERFAMILY MATCH")
                        else:
                            result_array = [line[0], line[6], scop_list[line[6]]]
                            results_list.append(result_array)
                    except:
                        result_array = [line[0], line[6], ""]
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
    result_dir = "/mnt/bioinf/archive0" \
                 "/eigen_thread/results/"
    start_parse_pattern = r"^\sNo\sHit"
    start_parse_re = re.compile(start_parse_pattern)
    end_parse_pattern = r"^No\s1"
    end_parse_re = re.compile(end_parse_pattern)
    length_pattern = r"^Match_columns\s(\d+)"
    length_re = re.compile(length_pattern)
    region_pattern = r"(\d+)-(\d+)"
    region_re = re.compile(region_pattern)
    for file in glob.glob(result_dir+"hhresults/*.hhr_modelled"):
        seq_length = float(0)
        print(file)
        results_list = []
        pdb = file[-18:-14]
        print(pdb)
        chain = file[-14]
        pdb_id = pdb+chain
        print(pdb_id)
        parse_ctl = False
        if pdb_id in bench_membership:
            print(result_dir+"processed_comparison_family/"+pdb_id+".hhtop")
            out = open(result_dir+"processed_comparison_family/"+pdb_id+".hhtop", "w+")
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
                            overlap_size = int(region_re_result.group(2))-int(region_re_result.group(1))
                        # print(scop_3_levels)

                        percentage = overlap_size/seq_length
                        if percentage < 0.4:
                            continue

                        try:
                            # print(scop_3_levels)
                            if scop_family == scop_class:
                                print("FAMILY MATCH")
                            # elif scop_3_levels == this_3_levels:
                                # print("SUPERFAMILY MATCH")
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
    result_dir = "/mnt/bioinf/archive0/" \
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
            print(result_dir+"processed_comparison_family/"+pdb_id+".genthtop")
            out = open(result_dir+"processed_comparison_family/"+pdb_id+".genthtop",
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
                    length = float(entries[8])
                    if overlap_size/length < 0.4:
                        continue
                    #
                    if scop_list[entries[9]] == scop_class:
                        print("FAMILY MATCH")
                    # elif scop_3_levels == this_3_levels:
                        # print("SUPERFAMILY MATCH")
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
        # break


parse_eigen(omit_from_results_set, scop_list, bench_membership)
# exit()
parse_hh(omit_from_results_set, scop_list, bench_membership)
parse_genth(omit_from_results_set_pdb, pdb_list, bench_membership)
