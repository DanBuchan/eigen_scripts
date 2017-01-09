import glob
from collections import defaultdict

pdb_list = []
for file in glob.glob('/mnt/bioinf/archive0/eigen_thread/eigenthreader/seq_files/*'):
    #print(file)
    pdb = file[-11:-7]
    pdb_list.append("d"+pdb+"a_")

count = 0
benchmarks = {}
with open("/mnt/bioinf/archive0/eigen_thread/scop_data/benchmark_domain_id.txt") as bench_ids:
    for line in bench_ids:
        line = line.rstrip()
        entries = line.split(",")
        benchmarks[entries[0]] = entries[1]

family_members = defaultdict(list)
superfamily_members = defaultdict(list)
fold_members = defaultdict(list)
with open("/mnt/bioinf/archive0/eigen_thread/scop_data/dir.des.scop.1.75.txt") as scop:
    for line in scop:
        line = line.rstrip()
        entries = line.split()
        if "px" in entries[1]:
            family = entries[2]
            superfamily = ".".join(entries[2].split(".")[:-1])
            fold = ".".join(superfamily.split(".")[:-1])
            pdbid = entries[3]
            family_members[family].append(pdbid)
            superfamily_members[superfamily].append(pdbid)
            fold_members[fold].append(pdbid)

# print(fold_members)
foldlib_members = []
for file in glob.glob("/scratch1/NOT_BACKED_UP/dbuchan/HHSearch_eigen/a3m/*"):
    domainid = file[-11:-4]
    # print(domainid)
    foldlib_members.append(domainid)

# print(foldlib_members)

def findFoldlibReps(membership, foldlib_members):
    for scopid in membership:
        intersect = set(membership[scopid]) & set(foldlib_members)
        # print(len(intersect))
        membership[scopid] = intersect
    return(membership)

family_members = findFoldlibReps(family_members, foldlib_members)
superfamily_members = findFoldlibReps(superfamily_members, foldlib_members)
fold_members = findFoldlibReps(fold_members, foldlib_members)

tp_family_removal = []
tp_superfamily_removal = []
for domain in benchmarks:
    domain_family = benchmarks[domain]
    domain_superfamily = ".".join(domain_family.split(".")[:-1])
    domain_fold = ".".join(domain_superfamily.split(".")[:-1])
    print(domain+","+domain_family+","+str(len(fold_members[domain_fold]))+"," +
          str(len(superfamily_members[domain_superfamily]))+"," +
          str(len(family_members[domain_family])) )
    if len(fold_members[domain_fold]) == len(family_members[domain_family]):
        tp_family_removal.append(domain)
    if len(fold_members[domain_fold]) == len(superfamily_members[domain_superfamily]):
        tp_superfamily_removal.append(domain)

print(tp_family_removal)
print(tp_superfamily_removal)

print(len(tp_family_removal))
print(len(tp_superfamily_removal))
