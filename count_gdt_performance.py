import glob
import pprint
import csv
import re
from subprocess import Popen, PIPE
import shlex
from statistics import mean

def skip_comments(iterable):
    for line in iterable:
        if not line.startswith('#'):
            yield line

result_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/" \
             "results/processed_comparison/"
eigen_models = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
               "eigen_thread/results/optimised_t20_c9/models/"
benchmark_models = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
                   "eigen_thread/eigenthreader/structures/"

tm_pattern = "TM-score\s+=\s+(.+)\s+\("
# gdt_pattern = "GDT-TS-score=\s+(.+?)\s+%"
tm_re = re.compile(tm_pattern)
# gdt_re = re.compile(gdt_pattern)

hh_tm_averages = [0, 0, 0, 0]
hh_gdt_averages = [0, 0, 0, 0]
gen_tm_averages = [0, 0, 0, 0]
gen_gdt_averages = [0, 0, 0, 0]


def average_eigen_scores():
    totes = 150
    eigen_tm_averages = [0, 0, 0, 0]
    eigen_gdt_averages = [0, 0, 0, 0]
    for file in glob.glob(result_dir+"*.eigentop"):
        pdb_id = file[-14:-9]
        print(file)
        with open(file) as scop_list_file:
            reader = csv.reader(skip_comments(scop_list_file), delimiter=',',
                                quotechar='"')
            tm_results = []
            gdt_results = []
            for line in reader:
                model = pdb_id+"_"+line[1]+".model.pdb"
                native_struct = pdb_id.upper()[0:4]+"_"+pdb_id.upper()[4:5]+".pdb"
                # print(benchmark_models+native_struct)
                # print(eigen_models+model)
                cmd = "/cs/research/bioinf/home1/green/dbuchan/bin/TMalign " + \
                    eigen_models+model+" " + \
                    benchmark_models+native_struct
                process = Popen(shlex.split(cmd), stdout=PIPE)
                (output, err) = process.communicate()
                if "There is no common residues in the input structures" in output:
                    print("Nope: "+pdb_id)
                    totes -= 1
                    continue
                exit_code = process.wait()
                # print(output.decode("utf-8") )
                print(cmd)
                tm_result = tm_re.search(output.decode("utf-8"))
                gdt_result = gdt_re.search(output.decode("utf-8"))

                tm_results.append(float(tm_result.group(1)))
                # gdt_results.append(float(gdt_result.group(1)))
                # print(tm_score)
                # > TMScore model native
                # break
            print(tm_results)
            print(gdt_results)

            eigen_tm_averages[0] += tm_results[0]
            eigen_gdt_averages[0] += gdt_results[0]
            eigen_tm_averages[1] += mean(tm_results[0:2])
            eigen_gdt_averages[1] += mean(gdt_results[0:2])
            eigen_tm_averages[2] += mean(tm_results[0:5])
            eigen_gdt_averages[2] += mean(gdt_results[0:5])
            eigen_tm_averages[3] += mean(tm_results)
            eigen_gdt_averages[3] += mean(gdt_results)


(eigen_tm_averages, eigen_gdt_averages) = average_eigen_scores()

print(eigen_tm_averages)
print(eigen_gdt_averages)
