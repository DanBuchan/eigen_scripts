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

def average_scores(result_dir, ending):

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

    tm_averages = [[], [], [], []]
    gdt_averages = [[], [], [], []]
    pdb_id = ''
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
        # print(pdb_id)
        # continue
        eprint(file)
        with open(file) as scop_list_file:
            reader = csv.reader(skip_comments(scop_list_file), delimiter=',',
                                quotechar='"')
            tm_results = []
            gdt_results = []
            for line in reader:
                if len(line[1]) == 0:
                    continue
                model = pdb_id+"_"+line[1]+".model.pdb"
                native_struct = pdb_id.upper()[0:4]+"_"+pdb_id.upper()[4:5]+".pdb"
                eprint(benchmark_models+native_struct)
                eprint(models_dir+model)
                # try:
                #     cmd = "/home/dbuchan/bin/TMalign " + \
                #         models_dir+model+" " + \
                #         benchmark_models+native_struct
                #     eprint(cmd)
                #     process = Popen(shlex.split(cmd), stdout=PIPE)
                #     (output, err) = process.communicate()
                #     #print(output.decode("utf-8"))
                #     exit_code = process.wait()
                #     tm_result = tm_re.search(output.decode("utf-8"))
                #     tm_results.append(float(tm_result.group(1)))
                # except:
                #     eprint("COULD NOT RUN TMalign")

                try:
                    cmd = "/home/dbuchan/bin/gdtlist " + \
                          benchmark_models+native_struct + \
                          " "+models_dir+model
                    eprint(cmd)
                    process = Popen(shlex.split(cmd), stdout=PIPE)
                    (output, err) = process.communicate()
                    exit_code = process.wait()
                    gdt_result = gdt_re.search(output.decode("utf-8"))
                    gdt_results.append(float(gdt_result.group(1)))
                    tm_result = tm_re.search(output.decode("utf-8"))
                    tm_results.append(float(tm_result.group(1)))
                except Exception as e:
                    eprint("COULD NOT RUN gdtlist")
                    eprint(str(e))
            #exit()

            # print(tm_results)
            # print(gdt_results)

            if len(tm_results) > 0:
                tm_averages[0].append(tm_results[0])
            if len(tm_results) > 1:
                tm_averages[1].append(max(tm_results[0:2]))
            if len(tm_results) > 4:
                tm_averages[2].append(max(tm_results[0:5]))
            if len(tm_results) > 9:
                tm_averages[3].append(max(tm_results))

            if len(gdt_results) > 0:
                gdt_averages[0].append(gdt_results[0])
            if len(gdt_results) > 1:
                gdt_averages[1].append(max(gdt_results[0:2]))
            if len(gdt_results) > 4:
                gdt_averages[2].append(max(gdt_results[0:5]))
            if len(gdt_results) > 9:
                gdt_averages[3].append(max(gdt_results))
        # break

    return(tm_averages, gdt_averages)

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

    tm_scores = defaultdict(list)
    gdt_scores = defaultdict(list)
    pdb_id = ''
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
            # try:
            #     cmd = "/home/dbuchan/bin/TMalign " + \
            #         models_dir+model+" " + \
            #         benchmark_models+native_struct
            #     eprint(cmd)
            #     process = Popen(shlex.split(cmd), stdout=PIPE)
            #     (output, err) = process.communicate()
            #     #print(output.decode("utf-8"))
            #     exit_code = process.wait()
            #     tm_result = tm_re.search(output.decode("utf-8"))
            #     tm_results.append(float(tm_result.group(1)))
            # except:
            #     eprint("COULD NOT RUN TMalign")

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
                    tm_scores[1].append(float(tm_result.group(1))
                    gdt_scores[1].append(float(gdt_result.group(1))
                if line_count <= 2:
                    tm_scores[2].append(float(tm_result.group(1))
                    gdt_scores[2].append(float(gdt_result.group(1))
                if line_count <= 5:
                    tm_scores[5].append(float(tm_result.group(1))
                    gdt_scores[5].append(float(gdt_result.group(1))
                if line_count <= 10:
                    tm_scores[10].append(float(tm_result.group(1))
                    gdt_scores[10].append(float(gdt_result.group(1))
            except Exception as e:
                eprint("COULD NOT RUN gdtlist")
                eprint(str(e))

    return(tm_scores, gdt_scores)
# (eigen_tm_averages, eigen_gdt_averages) = average_scores(result_dir,
#                                                          "*.eigentop")
# (gen_tm_averages, gen_gdt_averages) = average_scores(result_dir,
#                                                      "*.genthtop")
# (hh_tm_averages, hh_gdt_averages) = average_scores(result_dir,
#                                                    "*.hhtop")

(eigen_tm, eigen_gdt) = calculate_structure_scores(result_dir, "*.eigentop")
print(eigen_tm)
print(eigen_gdt)
# print(hh_tm_averages)
# print(hh_gdt_averages)
# print("level,eigen_max_tm,eigen_max_gdt,"
#       "gen_max_tm,gen_max_gdt,hh_max_tm,hh_max_gdt")
# print("t1,"+str(round(max(eigen_tm_averages[0]), 2))+"," +
#       str(round(max(eigen_gdt_averages[0]), 2))+"," +
#       str(round(max(gen_tm_averages[0]), 2))+"," +
#       str(round(max(gen_gdt_averages[0]), 2))+"," +
#       str(round(max(hh_tm_averages[0]), 2))+"," +
#       str(round(max(hh_gdt_averages[0]), 2)))
# print("t2,"+str(round(max(eigen_tm_averages[1]), 2))+"," +
#       str(round(max(eigen_gdt_averages[1]), 2))+"," +
#       str(round(max(gen_tm_averages[1]), 2))+"," +
#       str(round(max(gen_gdt_averages[1]), 2))+"," +
#       str(round(max(hh_tm_averages[1]), 2))+"," +
#       str(round(max(hh_gdt_averages[1]), 2)))
# print("t5,"+str(round(max(eigen_tm_averages[2]), 2))+"," +
#       str(round(max(eigen_gdt_averages[2]), 2))+"," +
#       str(round(max(gen_tm_averages[2]), 2))+"," +
#       str(round(max(gen_gdt_averages[2]), 2))+"," +
#       str(round(max(hh_tm_averages[2]), 2))+"," +
#       str(round(max(hh_gdt_averages[2]), 2)))
# print("t10,"+str(round(max(eigen_tm_averages[3]), 2))+"," +
#       str(round(max(eigen_gdt_averages[3]), 2))+"," +
#       str(round(max(gen_tm_averages[3]), 2))+"," +
#       str(round(max(gen_gdt_averages[3]), 2))+"," +
#       str(round(max(hh_tm_averages[3]), 2))+"," +
#       str(round(max(hh_gdt_averages[3]), 2)))
#
#
# print("level,eigen_median_tm,eigen_median_gdt,"
#       "gen_median_tm,gen_median_gdt,hh_median_tm,hh_median_gdt")
# print("t1,"+str(round(median(eigen_tm_averages[0]), 2))+"," +
#       str(round(median(eigen_gdt_averages[0]), 2))+"," +
#       str(round(median(gen_tm_averages[0]), 2))+"," +
#       str(round(median(gen_gdt_averages[0]), 2))+"," +
#       str(round(median(hh_tm_averages[0]), 2))+"," +
#       str(round(median(hh_gdt_averages[0]), 2)))
# print("t2,"+str(round(median(eigen_tm_averages[1]), 2))+"," +
#       str(round(median(eigen_gdt_averages[1]), 2))+"," +
#       str(round(median(gen_tm_averages[1]), 2))+"," +
#       str(round(median(gen_gdt_averages[1]), 2))+"," +
#       str(round(median(hh_tm_averages[1]), 2))+"," +
#       str(round(median(hh_gdt_averages[1]), 2)))
# print("t5,"+str(round(median(eigen_tm_averages[2]), 2))+"," +
#       str(round(median(eigen_gdt_averages[2]), 2))+"," +
#       str(round(median(gen_tm_averages[2]), 2))+"," +
#       str(round(median(gen_gdt_averages[2]), 2))+"," +
#       str(round(median(hh_tm_averages[2]), 2))+"," +
#       str(round(median(hh_gdt_averages[2]), 2)))
# print("t10,"+str(round(median(eigen_tm_averages[3]), 2))+"," +
#       str(round(median(eigen_gdt_averages[3]), 2))+"," +
#       str(round(median(gen_tm_averages[3]), 2))+"," +
#       str(round(median(gen_gdt_averages[3]), 2))+"," +
#       str(round(median(hh_tm_averages[3]), 2))+"," +
#       str(round(median(hh_gdt_averages[3]), 2)))
