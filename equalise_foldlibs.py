import glob
import os.path
import os

hhlib = '/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/HHSearch_hmm_complete/'
etlib = '/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/et_data/'

for file in glob.glob(hhlib+"*.hhm"):
    #print(file)
    domain = file[-11:-4]
    #print(domain)
    if not os.path.isfile(etlib+domain+".eig"):
        print("Missing "+domain)
        try:
            os.remove(file)
        except:
            print("hh file missing")
        try:
            os.remove(hhlib+domain+".a3m")
        except:
            print("a3m file missing")
