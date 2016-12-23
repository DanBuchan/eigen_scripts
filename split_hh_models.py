import glob
import re

start_pattern = r"^REMARK\s(\S+)\s(\w\.\d+\.\d+\.\d+)\s"
start_re = re.compile(start_pattern)

print_line = False
f2 = None

for file in glob.glob("/scratch0/NOT_BACKED_UP/dbuchan/hhresults/hhblits_corrected/*.pdb"):
    f1 = open(file, "r")
    pdb = file[-9:-4]
    print(pdb)
    for line in f1:
        start_result = start_re.match(line)
        if start_result:
            print_line = True
            #print(line)
            #print(pdb)
            f2 = open("/scratch0/NOT_BACKED_UP/dbuchan/hhresults/hhblits_corrected/models/"+pdb+"_"+start_result.group(1)+".model.pdb", "w")
        if print_line == True:
            f2.write(line)
        if "END" in line:
            f2.close()
            print_line=False
    #break
