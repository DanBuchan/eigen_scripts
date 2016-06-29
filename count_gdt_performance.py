
from __future__ import print_function
import glob
import pprint
import csv
import re
from subprocess import Popen, PIPE
import shlex
from statistics import mean
import sys


def skip_comments(iterable):
    for line in iterable:
        if not line.startswith('#'):
            yield line


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

result_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/" \
             "results/processed_comparison/"
eigen_models = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
               "eigen_thread/results/optimised_t20_c9/models/"
benchmark_models = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
                   "eigen_thread/eigenthreader/structures/"


def average_scores(result_dir, ending):

    if "eigentop" in ending:
        models_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
                     "eigen_thread/results/optimised_t20_c9/models/"
    if "genthtop" in ending:
        models_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
                     "eigen_thread/results/genthreader_results/models/"

    tm_pattern = "\nTM-score=\s+(.+?)\s+\("
    gdt_pattern = "\nGDT=\s+(.+)"
    tm_re = re.compile(tm_pattern)
    gdt_re = re.compile(gdt_pattern)

    tm_averages = [[], [], [], []]
    gdt_averages = [[], [], [], []]
    for file in glob.glob(result_dir+ending):
        pdb_id = file[-14:-9]
        eprint(file)
        with open(file) as scop_list_file:
            reader = csv.reader(skip_comments(scop_list_file), delimiter=',',
                                quotechar='"')
            tm_results = []
            gdt_results = []
            for line in reader:
                model = pdb_id+"_"+line[1]+".model.pdb"
                native_struct = pdb_id.upper()[0:4]+"_"+pdb_id.upper()[4:5]+".pdb"
                # print(benchmark_models+native_struct)
                # print(models_dir+model)
                try:
                    cmd = "/cs/research/bioinf/home1/green/dbuchan/bin/TMalign " + \
                        models_dir+model+" " + \
                        benchmark_models+native_struct
                    # print(cmd)
                    process = Popen(shlex.split(cmd), stdout=PIPE)
                    (output, err) = process.communicate()
                    #print(output.decode("utf-8"))
                    exit_code = process.wait()
                    tm_result = tm_re.search(output.decode("utf-8"))
                    tm_results.append(float(tm_result.group(1)))
                except:
                    eprint("COULD NOT RUN TMalign")

                try:
                    cmd = "/cs/research/bioinf/home1/green/dbuchan/bin/maxcluster64bit" + \
                        " -e "+models_dir+model + \
                        " -p "+benchmark_models+native_struct + \
                        " -in -gdt"
                    # print(cmd)
                    process = Popen(shlex.split(cmd), stdout=PIPE)
                    (output, err) = process.communicate()
                    exit_code = process.wait()
                    gdt_result = gdt_re.search(output.decode("utf-8"))
                    gdt_results.append(float(gdt_result.group(1)))
                except:
                    eprint("COULD NOT RUN maxcluster")

            # print(tm_results)
            # print(gdt_results)

            if len(tm_results) > 0:
                tm_averages[0].append(tm_results[0])
            if len(tm_results) > 1:
                tm_averages[1].append(mean(tm_results[0:2]))
            if len(tm_results) > 4:
                tm_averages[2].append(mean(tm_results[0:5]))
            if len(tm_results) > 9:
                tm_averages[3].append(mean(tm_results))

            if len(gdt_results) > 0:
                gdt_averages[0].append(gdt_results[0])
            if len(gdt_results) > 1:
                gdt_averages[1].append(mean(gdt_results[0:2]))
            if len(gdt_results) > 4:
                gdt_averages[2].append(mean(gdt_results[0:5]))
            if len(gdt_results) > 9:
                gdt_averages[3].append(mean(gdt_results))

    return(tm_averages, gdt_averages)

(eigen_tm_averages, eigen_gdt_averages) = average_scores(result_dir,
                                                         "*.eigentop")
(gen_tm_averages, gen_gdt_averages) = average_scores(result_dir,
                                                     "*.genthtop")

tm_results = []
gdt_results = []
tm_pattern = "\nTM-score=\s+(.+?)\s+\("
gdt_pattern = "\nGDT=\s+(.+)"
tm_re = re.compile(tm_pattern)
gdt_re = re.compile(gdt_pattern)
for file in glob.glob("/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/results/hhresults/models/"+"*.pdb"):
    pdb_id = file[-9:-4]
    native_struct = pdb_id.upper()[0:4]+"_"+pdb_id.upper()[4:5]+".pdb"
    try:
        cmd = "/cs/research/bioinf/home1/green/dbuchan/bin/TMalign " + \
            file+" " + \
            benchmark_models+native_struct
        # print(cmd)
        process = Popen(shlex.split(cmd), stdout=PIPE)
        (output, err) = process.communicate()
        #print(output.decode("utf-8"))
        exit_code = process.wait()
        tm_result = tm_re.search(output.decode("utf-8"))
        tm_results.append(float(tm_result.group(1)))
    except:
        eprint("COULD NOT RUN TMalign")

    try:
        cmd = "/cs/research/bioinf/home1/green/dbuchan/bin/maxcluster64bit" + \
            " -e "+file + \
            " -p "+benchmark_models+native_struct + \
            " -in -gdt"
        # print(cmd)
        process = Popen(shlex.split(cmd), stdout=PIPE)
        (output, err) = process.communicate()
        exit_code = process.wait()
        gdt_result = gdt_re.search(output.decode("utf-8"))
        gdt_results.append(float(gdt_result.group(1)))
    except:
        eprint("COULD NOT RUN maxcluster")
    #break

print("level,eigen_average_tm,eigen_average_gdt,"
      "gen_average_tm,gen_average_gdt,hh_average_tm,hh_average_gdt")
print("t1,"+str(round(mean(eigen_tm_averages[0]), 2))+"," +
      str(round(mean(eigen_gdt_averages[0]), 2))+"," +
      str(round(mean(gen_tm_averages[0]), 2))+"," +
      str(round(mean(gen_gdt_averages[0]), 2))+",-,-")
print("t2,"+str(round(mean(eigen_tm_averages[1]), 2))+"," +
      str(round(mean(eigen_gdt_averages[1]), 2))+"," +
      str(round(mean(gen_tm_averages[1]), 2))+"," +
      str(round(mean(gen_gdt_averages[1]), 2))+",-,-")
print("t5,"+str(round(mean(eigen_tm_averages[2]), 2))+"," +
      str(round(mean(eigen_gdt_averages[2]), 2))+"," +
      str(round(mean(gen_tm_averages[2]), 2))+"," +
      str(round(mean(gen_gdt_averages[2]), 2))+",-,-")
print("t10,"+str(round(mean(eigen_tm_averages[3]), 2))+"," +
      str(round(mean(eigen_gdt_averages[3]), 2))+"," +
      str(round(mean(gen_tm_averages[3]), 2))+"," +
      str(round(mean(gen_gdt_averages[3]), 2))+"," +
      str(round(mean(tm_results), 2))+","+str(round(mean(gdt_results), 2)))
