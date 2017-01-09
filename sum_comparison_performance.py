# SCOP classification
#
# CLASS.fold.superfamilies.families

import glob
import pprint
import re
import math
from collections import defaultdict
pp = pprint.PrettyPrinter(indent=4)

family_removal_set = ['1i4jA', '1j3aA', '1chdA', '1dmgA', '1mk0A', '1wjxA',
                      '1jkxA', '1lm4A', '1g2rA', '1k6kA', '1fcyA', '3dqgA',
                      '1xkrA', '1dixA', '1nrvA', '1hxnA', '1dbxA', '1whiA',
                      '1im5A', '1kq6A']
superfamily_removal_set = ['1ej0A', '1ne2A', '1i4jA', '1gbsA', '5ptpA',
                           '1xdzA', '1a3aA', '1k7jA', '1j3aA', '1fk5A',
                           '1aoeA', '1behA', '1kqrA', '1i58A', '1chdA',
                           '3borA', '1iwdA', '1dmgA', '1mk0A', '1npsA',
                           '1gzcA', '1htwA', '1wjxA', '1ktgA', '1i1nA',
                           '1jkxA', '1g9oA', '2hs1A', '1jbkA', '1ql0A',
                           '1lm4A', '1g2rA', '1atzA', '1f6bA', '1gz2A',
                           '1i71A', '1k6kA', '1bebA', '1hdoA', '1fcyA',
                           '3dqgA', '1fqtA', '1dsxA', '1jyhA', '1h2eA',
                           '1cxyA', '1xkrA', '1vhuA', '1d4oA', '1dixA',
                           '1nrvA', '1vfyA', '1lo7A', '1hxnA', '1ckeA',
                           '1cjwA', '1c44A', '1wkcA', '1qf9A', '1dbxA',
                           '1tzvA', '1h0pA', '1c52A', '1aapA', '1ctfA',
                           '1gmxA', '1bsgA', '1ek0A', '1ihzA', '1mugA',
                           '1whiA', '1im5A', '1kq6A', '1tqhA']

def process_results(t1_results, t2_results, t5_results, t10_results,
                    file_ending, result_dir, removal_set):
    pdb_pattern = r"#\sPDB\sID:\s(.+)"
    scop_pattern = r"#\sSCOP\sFAMILY:\s(.+)"
    result_pattern = r".+,.+,.+"
    pdb_re = re.compile(pdb_pattern)
    scop_re = re.compile(scop_pattern)
    result_re = re.compile(result_pattern)
    this_pdb = ''
    this_scop_class = ''
    scop_class = ''
    fold = ''
    superf = ''
    family = ''

    param = file_ending
    t1_results[param]["class"] = 0
    t2_results[param]["class"] = 0
    t5_results[param]["class"] = 0
    t10_results[param]["class"] = 0
    t1_results[param]["fold"] = 0
    t2_results[param]["fold"] = 0
    t5_results[param]["fold"] = 0
    t10_results[param]["fold"] = 0
    t1_results[param]["superf"] = 0
    t2_results[param]["superf"] = 0
    t5_results[param]["superf"] = 0
    t10_results[param]["superf"] = 0
    t1_results[param]["family"] = 0
    t2_results[param]["family"] = 0
    t5_results[param]["family"] = 0
    t10_results[param]["family"] = 0

    # print(result_dir+"*."+file_ending)
    for file in glob.glob(result_dir+"*."+file_ending):
        #print(file)
        pdb = file.rsplit('/', 1)[-1][0:5]
        #print(pdb)
        if pdb in removal_set:
            continue
        counts = {}
        # print(file[-9:])
        result_count = 0
        f1 = open(file, "r")
        class_found = 100
        fold_found = 100
        superf_found = 100
        family_found = 100
        for line in f1:
            line = line.rstrip()
            line = line.rstrip(" (")
            # print(line)
            pdb_result = pdb_re.match(line)
            scop_result = scop_re.match(line)
            line_result = result_re.match(line)
            if(pdb_result):
                this_pdb = pdb_result.group(1)
            elif(scop_result):
                this_scop_class = scop_result.group(1)
                # print(this_scop_class)
                scop_class, fold, this_superf, this_family = this_scop_class.split(".")
            elif(line_result):
                result_count += 1
                entries = line.split(",")
                entries = entries[2:]
                # print(entries)
                for entry in entries:
                    #print(entry)
                    this_class, this_fold, this_superf, this_family = entry.split(".")
                    if this_class == scop_class and result_count < class_found:
                        class_found = result_count
                    if this_class == scop_class and this_fold == fold and \
                       result_count < fold_found:
                        fold_found = result_count
                    if this_class == scop_class and this_fold == fold and \
                       this_superf == superf and result_count < superf_found:
                        superf_found = result_count
                    if this_class == scop_class and this_fold == fold and \
                       this_superf == superf and this_family == family and \
                       result_count < family_found:
                        family_found = result_count

        if class_found == 1:
            t1_results[param]["class"] += 1
            t2_results[param]["class"] += 1
            t5_results[param]["class"] += 1
            t10_results[param]["class"] += 1
        if class_found == 2:
            t2_results[param]["class"] += 1
            t5_results[param]["class"] += 1
            t10_results[param]["class"] += 1
        if class_found <= 5 and class_found > 2:
            t5_results[param]["class"] += 1
            t10_results[param]["class"] += 1
        if class_found <= 10 and class_found > 5:
            t10_results[param]["class"] += 1

        if fold_found == 1:
            t1_results[param]["fold"] += 1
            t2_results[param]["fold"] += 1
            t5_results[param]["fold"] += 1
            t10_results[param]["fold"] += 1
        if fold_found == 2:
            t2_results[param]["fold"] += 1
            t5_results[param]["fold"] += 1
            t10_results[param]["fold"] += 1
        if fold_found <= 5 and fold_found > 2:
            t5_results[param]["fold"] += 1
            t10_results[param]["fold"] += 1
        if fold_found <= 10 and fold_found > 5:
            t10_results[param]["fold"] += 1

        if superf_found == 1:
            t1_results[param]["superf"] += 1
            t2_results[param]["superf"] += 1
            t5_results[param]["superf"] += 1
            t10_results[param]["superf"] += 1
        if superf_found == 2:
            t2_results[param]["superf"] += 1
            t5_results[param]["superf"] += 1
            t10_results[param]["superf"] += 1
        if superf_found <= 5 and superf_found > 2:
            t5_results[param]["superf"] += 1
            t10_results[param]["superf"] += 1
        if superf_found <= 10 and superf_found > 5:
            t10_results[param]["superf"] += 1

        #print(family_found)
        if family_found == 1:
            t1_results[param]["family"] += 1
            t2_results[param]["family"] += 1
            t5_results[param]["family"] += 1
            t10_results[param]["family"] += 1
        if family_found == 2:
            t2_results[param]["family"] += 1
            t5_results[param]["family"] += 1
            t10_results[param]["family"] += 1
        if family_found <= 5 and family_found > 2:
            t5_results[param]["family"] += 1
            t10_results[param]["family"] += 1
        if family_found <= 10 and family_found > 5:
            t10_results[param]["family"] += 1

        # break

    return(t1_results, t2_results, t5_results, t10_results)

def print_counts(top, data, removal_set):
    output_line = ''
    denominator = 150 - len(removal_set)
    for vect_count in sorted(data.items()):
        output_line += top+","+str(vect_count[0])+","
        output_line += str(round(float(vect_count[1]["class"])/denominator, 3))+","
        output_line += str(round(float(vect_count[1]["fold"])/denominator, 3))+","
        output_line += str(round(float(vect_count[1]["superf"])/denominator, 3))+","
        output_line += str(round(float(vect_count[1]["family"])/denominator, 3))
        output_line += "\n"
    print(output_line.rstrip())


def process_data_set(result_dir, removal_set):

    t1_results = defaultdict(dict)
    t2_results = defaultdict(dict)
    t5_results = defaultdict(dict)
    t10_results = defaultdict(dict)
    t1_results, t2_results, t5_results, t10_results = process_results(t1_results,
                                                                      t2_results,
                                                                      t5_results,
                                                                      t10_results,
                                                                      "eigentop",
                                                                      result_dir,
                                                                      removal_set,
                                                                      )
    t1_results, t2_results, t5_results, t10_results = process_results(t1_results,
                                                                      t2_results,
                                                                      t5_results,
                                                                      t10_results,
                                                                      "hhtop",
                                                                      result_dir,
                                                                      removal_set,
                                                                      )
    t1_results, t2_results, t5_results, t10_results = process_results(t1_results,
                                                                      t2_results,
                                                                      t5_results,
                                                                      t10_results,
                                                                      "genthtop",
                                                                      result_dir,
                                                                      removal_set,
                                                                      )
    # pp.pprint(t1_results)

    print("top,vectors,class,fold,superf,family")
    print_counts("1", t1_results, removal_set)
    print_counts("2", t2_results, removal_set)
    print_counts("5", t5_results, removal_set)
    print_counts("10", t10_results, removal_set)

process_data_set("/mnt/bioinf/archive0/eigen_thread/results/processed_comparison_family/", family_removal_set)
#process_data_set("/mnt/bioinf/archive0/eigen_thread/results/processed_comparison_superfamily/", superfamily_removal_set)
