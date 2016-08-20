import glob
import csv
import pprint
from random import shuffle

pp = pprint.PrettyPrinter(indent=4)


def file_printer(list_size, file, contact_list):
    f = open(file, 'w')
    for i in range(0, int(list_size)+1):
        f.write(contact_list[i][0]+" "+contact_list[i][1]+" " +
                contact_list[i][2]+" "+contact_list[i][3]+" " +
                contact_list[i][4]+"\n")

# build a look up of seq lengths
length_lookup = dict()
seq_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/" \
          "eigenthreader/seq_files/*.fasta"
for file in glob.glob(seq_dir):
    with open(file) as fasta_file:
        pdb = file[-11:-6]
        # print(pdb)
        for line in fasta_file:
            if ">" in line:
                continue
            # print(len(line))
            length_lookup[pdb] = len(line)

# print(length_lookup)
# read a file

contact_dir = "/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/" \
              "eigenthreader/confiles/*.con"

for file in glob.glob(contact_dir):
    pdb = file[-9:-4]
    # print(pdb)
    contact_list = []
    contact_sorted = []
    contact_rand = []
    with open(file) as con_file:
        reader = csv.reader(con_file, delimiter=' ',
                            quotechar='"')
        pdb_length = length_lookup[pdb]
        for line in reader:
            if abs(int(line[0])-int(line[1])) > 21:
                contact_list.append(line)
        contact_sorted = sorted(contact_list, key=lambda x: float(x[4]),
                                reverse=True)

    file_printer(pdb_length, "/cs/research/bioinf/home1/green/dbuchan/"
                             "archive0/eigen_thread/contact_subsets/top/"
                             "L/"+pdb+"_L.con", contact_sorted)
    file_printer(pdb_length/2, "/cs/research/bioinf/home1/green/dbuchan/"
                               "archive0/eigen_thread/contact_subsets/top/"
                               "L2/"+pdb+"_L2.con", contact_sorted)
    file_printer(pdb_length/5, "/cs/research/bioinf/home1/green/dbuchan/"
                               "archive0/eigen_thread/contact_subsets/top/"
                               "L5/"+pdb+"_L5.con", contact_sorted)
    file_printer(pdb_length/10, "/cs/research/bioinf/home1/green/dbuchan/"
                                "archive0/eigen_thread/contact_subsets/top/"
                                "L10/"+pdb+"_L10.con", contact_sorted)
    shuffle(contact_sorted)
    file_printer(pdb_length, "/cs/research/bioinf/home1/green/dbuchan/"
                             "archive0/eigen_thread/contact_subsets/random/"
                             "L/"+pdb+"_L.con", contact_sorted)
    shuffle(contact_sorted)
    file_printer(pdb_length/2, "/cs/research/bioinf/home1/green/dbuchan/"
                               "archive0/eigen_thread/contact_subsets/random/"
                               "L2/"+pdb+"_L2.con", contact_sorted)
    shuffle(contact_sorted)
    file_printer(pdb_length/5, "/cs/research/bioinf/home1/green/dbuchan/"
                               "archive0/eigen_thread/contact_subsets/random/"
                               "L5/"+pdb+"_L5.con", contact_sorted)
    shuffle(contact_sorted)
    file_printer(pdb_length/10, "/cs/research/bioinf/home1/green/dbuchan/"
                                "archive0/eigen_thread/contact_subsets/random/"
                                "L10/"+pdb+"_L10.con", contact_sorted)
    # break
