import glob
import pprint
import csv
import re


def skip_comments(iterable):
    for line in iterable:
        if not line.startswith('#'):
            yield line

result_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/" \
             "results/processed_comparison/"
eigen_models = "/cs/research/bioinf/home1/green/dbuchan/archive0/" \
               "eigen_thread/results/optimised_t20_c9/models/"

for file in glob.glob(result_dir+"*.eigentop"):
    pdb_id = file[-14:-9]
    print(file)
    with open(file) as scop_list_file:
        reader = csv.reader(skip_comments(scop_list_file), delimiter=',',
                            quotechar='"')
        for line in reader:
            model = pdb_id+"_"+line[1]+".model.pdb"
            print(model)
    break
