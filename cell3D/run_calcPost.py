#!/usr/bin/python 

import shutil
import os
import sys

def copyParameters(dirName):
   dirName = dirName+"/"
   shutil.copy(dirName+"customConstants.h",os.getcwd())
   shutil.copy(dirName+"parameters.h",os.getcwd())

def makePostProcess():
   os.system("make clean")
   os.system("make -j 8 calcPost")

def runPostProcess(dirName,name,start,stop):
   currDir = os.getcwd()
   os.chdir(dirName)
   
   filename = "pbs-script"
   file = open(filename,"w")
   file.write("#PBS -N calcPost\n")  
   file.write("#PBS -l nodes=1:ppn=1\n")
   file.write(open("pbs-script-base","r").read())
   file.write("./calcPost "+name+" "+str(start)+" "+str(stop) +" \n")
   file.close()

   os.system("qsub pbs-script")

   os.chdir(currDir)

def copyToDirectory(dirName):
   dirName = dirName+"/"
   shutil.copy("calcPost",dirName)
   shutil.copy("pbs-script-base",dirName)

dir = "/home/mrobinson/data/cell3D/120801174033variable_liqdemTS_testAll"
subdirs3 = [["d1_v100_dry2cpu8",99,274]]
subdirs = [["d1_v100_dry2cpu8",99,274], ["d1_v400_dry2cpu8",99,123], ["d2_v100_dry5cpu8",99,163], ["d2_v400_dry5cpu8",99,256], ["d5_v100_dry14cpu8",99,166], ["d5_v400_dry14cpu8",99,327]]
subdirs2 = [["d1_v50_dry2cpu8",99,481],["d1_v200_dry2cpu8",99,113],["d1_v600_dry2cpu8",100,107]] 

subdirs = subdirs + subdirs2
#subdirs = subdirs3
#subdirs = [["d1_v100_dry2cpu8",99,274]]
name = "results"

for subdir in subdirs:
   dirName = dir+"/"+subdir[0]
   copyParameters(dirName)
   makePostProcess()
   copyToDirectory(dirName)
   runPostProcess(dirName,name,subdir[1],subdir[2])
