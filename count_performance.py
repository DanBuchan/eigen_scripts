# SCOP classification
#
# CLASS.fold.superfamilies.families

import glob
import pprint
import re
import math
from collections import defaultdict
pp = pprint.PrettyPrinter(indent=4)

result_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/results/random"
pdb_pattern = r"#\sPDB\sID:\s(.+)"
scop_pattern = r"#\sSCOP\sFAMILY:\s(.+)"
result_pattern = r".+,.+,(.+)"
pdb_re = re.compile(pdb_pattern)
scop_re = re.compile(scop_pattern)
result_re = re.compile(result_pattern)
this_pdb = ''
this_scop_class = ''
scop_class = ''
fold = ''
superf = ''
family = ''
t1_results = defaultdict(dict)
t2_results = defaultdict(dict)
t5_results = defaultdict(dict)
t10_results = defaultdict(dict)

for directory in glob.glob(result_dir+"/L*"):
    param = directory[-3:]
    param = int(param.lstrip("/L"))
    # print(param)
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
    # print(param)
    for file in glob.glob(directory+"/*.top"):
        counts = {}
        # print(file[-9:])
        result_count = 0
        f1 = open(file, "r")
        class_found = 100
        fold_found = 100
        superf_found = 100
        family_found = 100
        for line in f1:
            pdb_result = pdb_re.match(line)
            scop_result = scop_re.match(line)
            line_result = result_re.match(line)
            if(pdb_result):
                this_pdb = pdb_result.group(1)
            elif(scop_result):
                this_scop_class = scop_result.group(1)
                scop_class, fold, superf, family = this_scop_class.split(".")
            elif(line_result):
                result_count += 1
                this_class, this_fold, this_superf, this_family = line_result.group(1).split(".")
                if this_class == scop_class and result_count < class_found:
                    class_found = result_count
                if this_class == scop_class and this_fold == fold and result_count < fold_found:
                    fold_found = result_count
                if this_class == scop_class and this_fold == fold and  this_superf == superf and result_count < superf_found:
                    superf_found = result_count
                if this_class == scop_class and this_fold == fold and  this_superf == superf and  this_family == family and result_count < family_found:
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
    # break

#print(t10_results)

def print_counts(top, data):
    output_line = ''
    for vect_count in sorted(data.items()):
            output_line += top+","+str(vect_count[0])+","
            output_line += str(round(vect_count[1]["class"]/150, 3))+","
            output_line += str(round(vect_count[1]["fold"]/150, 3))+","
            output_line += str(round(vect_count[1]["superf"]/150, 3))+","
            output_line += str(round(vect_count[1]["family"]/150, 3))
            output_line += "\n"
    print(output_line.rstrip())
print("top,vectors,class,fold,superf,family")
print_counts("1", t1_results)
print_counts("2", t2_results)
print_counts("5", t5_results)
print_counts("10", t10_results)
