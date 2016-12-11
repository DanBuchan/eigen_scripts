import glob
import csv
from collections import defaultdict

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

# scop_list we now have a list of scop IDs (dxxxxx) to scop classes (x.x.x.x)

bench_prots = {}
bench_membership = {}
with open("/mnt/bioinf/archive0/eigen_thread/scop_data/benchmark_family_members.txt") as benchlist:
    reader = csv.reader(benchlist, delimiter=',', quotechar='"')
    for row in reader:
        bench_prots[row[0]] = row[2:len(row)]
        bench_membership[row[0]] = row[1]

# with open("/mnt/bioinf/archive0/eigen_thread/scop_data/benchmark_superfamily_members.txt") as benchlist:
#     reader = csv.reader(benchlist, delimiter=',', quotechar='"')
#     for row in reader:
#         bench_prots[row[0]]+row[2:len(row)]
#         bench_membership[row[0]] = row[1]

for prot in bench_prots:
    bench_prots[prot] = list(set(bench_prots[prot]))

print("Bench members"+str(len(bench_prots["1abaA"])))
# bench_prots a list of the benchmark proteins with all their superfamily and family members
# benchmembershop a list of benchmark pdb ID to scop class (x.x.x.x)

# print(bench_prots["1abaA"])
fold_lib_remaining = defaultdict(list)
# count_folds = 0
for pdb in bench_prots:
    if "1abaA" not in pdb:
        continue
    print(pdb)
    for file in glob.glob("/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/et_foldlib/*.eig"):
        scopid = file[-11:-4]
        if scopid not in bench_prots[pdb]:
            if scopid is pdb:
                continue
            fold_lib_remaining[pdb].append(scopid)
            # if "c.47." in scop_list[scopid]:
            #     count_folds += 1

# print(count_folds)
print(len(fold_lib_remaining["1abaA"]))
# fold_lib_remaining, for each benchmark pdb id we have a list of all the
# folds in the fold_lib that WERE NOT in the superfamily/family set.
# now count up howmany of the proteins have the same fold_lib
fold_count = defaultdict(lambda:0)
for pdb in fold_lib_remaining:
    pdb_fold = ".".join(bench_membership[pdb].split(".")[0:2])
    for scop_id in fold_lib_remaining[pdb]:
        member_fold = ".".join(scop_list[scop_id].split(".")[0:2])
        print(pdb+" "+pdb_fold+" : "+scop_id+" "+member_fold)
        if pdb_fold in member_fold:
             fold_count[pdb] += 1
        # # print(pdb_fold+" : "+member_fold)

print(fold_count)
