from __future__ import print_function
import glob
import pprint
import csv
import re
from subprocess import Popen, PIPE
import shlex
from statistics import median
import sys
from collections import defaultdict


def add_values(key, max_scores, median_scores, values_list):
    try:
        max_scores[key].append(round(max(values_list[key]), 3))
    except Exception as e:
        eprint("Could not add max values to results")
        eprint(str(e))
    try:
        median_scores[key].append(round(median(values_list[key]), 3))
    except Exception as e:
        eprint("Could not add median values to results")
        eprint(str(e))
    return(max_scores, median_scores)

def skip_comments(iterable):
    for line in iterable:
        if not line.startswith('#'):
            yield line


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

result_dir = "/mnt/bioinf/archive0/eigen_thread/" \
             "results/processed_comparison/"
eigen_models = "/mnt/bioinf/archive0/" \
               "eigen_thread/results/optimised_t20_c9/models/"
benchmark_models = "/mnt/bioinf/archive0/" \
                   "eigen_thread/eigenthreader/structures/"


# def getFastaLength(file):
#     with open(file) as fasta_file:
#         for line in fasta_file:
#             if ">" in line:
#                 continue
#             else:
#                 return(len(line.rstrip()))

def read_missing(file):
    missing_things = []
    with open(file) as missing_file:
        for line in missing_file:
            line = line.rstrip()
            missing_things.append(line)
    return(missing_things)

def calculate_structure_scores(result_dir, ending):

    if "eigentop" in ending:
        models_dir = "/mnt/bioinf/archive0/" \
                     "eigen_thread/results/optimised_t20_c9/models/"
    if "genthtop" in ending:
        models_dir = "/mnt/bioinf/archive0/" \
                     "eigen_thread/results/genthreader_results/models/"
    if "hhtop" in ending:
        models_dir = "/mnt/bioinf/archive0/" \
                     "eigen_thread/results/hhresults/models/"

    missing_things = read_missing("/home/dbuchan/Code/eigen_scripts/"
                                  "missing_hh_models.txt")

    tm_pattern = "TM-Score\s=\s(.+?)\sRMSD"
    gdt_pattern = "GDT-TS\s=\s(.+)\sTM"
    tm_re = re.compile(tm_pattern)
    gdt_re = re.compile(gdt_pattern)

    tm_max_scores = defaultdict(list)
    gdt_max_scores = defaultdict(list)
    tm_median_scores = defaultdict(list)
    gdt_median_scores = defaultdict(list)

    pdb_id = ''
    tm_t1_scores = {}
    gdt_t1_scores = {}

    for file in glob.glob(result_dir+ending):
        if "eigentop" in ending or "genthtop" in ending:
            pdb_id = file[-14:-9]
        else:
            pdb_id = file[-11:-6]
        # fasta_file = "/mnt/bioinf/archive0/" \
        #              "eigen_thread/eigenthreader/seq_files/"+pdb_id+".fasta"
        #fasta_length = getFastaLength(fasta_file)
        #print(pdb_id)
        if pdb_id in missing_things:
            print("MISSING")
            continue

        eprint(file)
        t1_tm = defaultdict(list)
        t2_tm = defaultdict(list)
        t5_tm = defaultdict(list)
        t10_tm = defaultdict(list)
        t1_gdt = defaultdict(list)
        t2_gdt = defaultdict(list)
        t5_gdt = defaultdict(list)
        t10_gdt = defaultdict(list)
        with open(file) as scop_list_file:
            reader = csv.reader(skip_comments(scop_list_file), delimiter=',',
                                quotechar='"')
            line_count = 0
            for line in reader:
                if len(line[1]) == 0:
                    continue
                model = pdb_id+"_"+line[1]+".model.pdb"
                native_struct = pdb_id.upper()[0:4]+"_"+pdb_id.upper()[4:5]+".pdb"
                eprint(benchmark_models+native_struct)
                eprint(models_dir+model)

                try:
                    cmd = "/home/dbuchan/bin/gdtlist " + \
                          benchmark_models+native_struct + \
                          " "+models_dir+model
                    eprint(cmd)
                    process = Popen(shlex.split(cmd), stdout=PIPE)
                    (output, err) = process.communicate()
                    exit_code = process.wait()
                    line_count+=1
                    gdt_result = gdt_re.search(output.decode("utf-8"))
                    tm_result = tm_re.search(output.decode("utf-8"))
                    if line_count == 1:
                        t1_tm[1].append(float(tm_result.group(1)))
                        t1_gdt[1].append(float(gdt_result.group(1)))
                        tm_t1_scores[pdb_id]=float(tm_result.group(1))
                    if line_count <= 2:
                        t2_tm[2].append(float(tm_result.group(1)))
                        t2_gdt[2].append(float(gdt_result.group(1)))
                    if line_count <= 5:
                        t5_tm[5].append(float(tm_result.group(1)))
                        t5_gdt[5].append(float(gdt_result.group(1)))
                    if line_count <= 10:
                        t10_tm[10].append(float(tm_result.group(1)))
                        t10_gdt[10].append(float(gdt_result.group(1)))
                except Exception as e:
                    eprint("COULD NOT RUN gdtlist: ")
                    eprint(str(e))

        (tm_max_scores, tm_median_scores) = add_values(1, tm_max_scores,
                                                       tm_median_scores, t1_tm)
        (tm_max_scores, tm_median_scores) = add_values(2, tm_max_scores,
                                                       tm_median_scores, t2_tm)
        (tm_max_scores, tm_median_scores) = add_values(5, tm_max_scores,
                                                       tm_median_scores, t5_tm)
        (tm_max_scores, tm_median_scores) = add_values(10, tm_max_scores,
                                                       tm_median_scores, t10_tm)
        (gdt_max_scores, gdt_median_scores) = add_values(1, gdt_max_scores,
                                                         gdt_median_scores,
                                                         t1_gdt)
        (gdt_max_scores, gdt_median_scores) = add_values(2, gdt_max_scores,
                                                         gdt_median_scores,
                                                         t2_gdt)
        (gdt_max_scores, gdt_median_scores) = add_values(5, gdt_max_scores,
                                                         gdt_median_scores,
                                                         t5_gdt)
        (gdt_max_scores, gdt_median_scores) = add_values(10, gdt_max_scores,
                                                         gdt_median_scores,
                                                         t10_gdt)

    return(tm_max_scores, tm_median_scores, gdt_max_scores, gdt_median_scores,
           tm_t1_scores)

(eigen_max_tm,
 eigen_median_tm,
 eigen_max_gdt,
 eigen_median_gdt,
 eigen_tm_t1_scores) = calculate_structure_scores(result_dir, "*.eigentop")
(gen_max_tm,
 gen_median_tm,
 gen_max_gdt,
 gen_median_gdt,
 gen_tm_t1_scores) = calculate_structure_scores(result_dir, "*.genthtop")
(hh_max_tm,
 hh_median_tm,
 hh_max_gdt,
 hh_median_gdt,
 hh_tm_t1_scores) = calculate_structure_scores(result_dir, "*.hhtop")

# print(eigen_max_tm)
# print(eigen_median_tm)
# print(eigen_max_gdt)
# print(eigen_median_gdt)

print("tm_level,eigen_max,eigen_median,gen_max,gen_median,hh_max,hh_median")
for level in sorted(eigen_max_tm):
    info_str = str(level)+","+str(round(median(eigen_max_tm[level]), 2))
    info_str += ","+str(round(median(eigen_median_tm[level]), 2))
    info_str += ","+str(round(median(gen_max_tm[level]), 2))
    info_str += ","+str(round(median(gen_median_tm[level]), 2))
    info_str += ","+str(round(median(hh_max_tm[level]), 2))
    info_str += ","+str(round(median(hh_median_tm[level]), 2))
    print(info_str)
print("gdt_level,eigen_max,eigen_median,gen_max,gen_median,hh_max,hh_median")
for level in sorted(eigen_max_gdt):
    info_str = str(level)+","+str(round(median(eigen_max_gdt[level]), 2))
    info_str += ","+str(round(median(eigen_median_gdt[level]), 2))
    info_str += ","+str(round(median(gen_max_gdt[level]), 2))
    info_str += ","+str(round(median(gen_median_gdt[level]), 2))
    info_str += ","+str(round(median(hh_max_gdt[level]), 2))
    info_str += ","+str(round(median(hh_median_gdt[level]), 2))
    print(info_str)


for pdb in hh_tm_t1_scores:
    print(pdb+",hh,"+str(hh_tm_t1_scores[pdb]))
    try:
        print(pdb+",eigen,"+str(eigen_tm_t1_scores[pdb]))
    except:
        print(pdb+",eigen,-")
