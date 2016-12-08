import glob
import csv

bench_prots = {}
with open("/mnt/bioinf/archive0/eigen_thread/scop_data/benchmark_superfamily_members.txt") as benchlist:
    reader = csv.reader(benchlist, delimiter=',', quotechar='"')
    for row in reader:
        #print(row[0])
        #nprint(row[2:len(row)])
        bench_prots[row[0]] = row[2:len(row)]

#print(bench_prots)

for file in glob.glob("/scratch0/NOT_BACKED_UP/dbuchan/eigen_benchmark/et_foldlib/*.eig"):
    #print(file)
    scopid= file[-11:-4]
    for pdb in bench_prots:
        if scopid in bench_prots[pdb]:
            bench_prots[pdb] = [x for x in bench_prots[pdb] if x != scopid]

print(bench_prots)

for pdb in bench_prots:
    print(pdb+" : "+str(len()))
