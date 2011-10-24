#!/usr/bin/python 
import sys
import os
import shutil,glob
import datetime

def writeParametersMulti(p,res,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define MANY_PARTICLES\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.close()


def writeParametersMultiGrid(p,res,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define MANY_PARTICLES\n")
   file.write("#define GRID_OF_DEM\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.close()


def myRemove(name):
   if os.path.exists(name):
      os.remove(name)

def compileProgram():
   os.system("make clean")
   os.system("make -j 8")

def copyPattern(pattern,dest_dir):
   files = glob.iglob(pattern)
   for file in files:
      if os.path.isfile(file):
         shutil.copy2(file, dest_dir)

def copyToDirectory(dirName):
   d = os.path.dirname(dirName);
   if not os.path.exists(dirName):
      os.makedirs(dirName)
   dirName = dirName+"/"
   shutil.copy("run",dirName)
   shutil.copy("setup",dirName)
   shutil.copy("createPVD",dirName)
   shutil.copy("customConstants.h",dirName)
   shutil.copy("parameters.h",dirName)
   copyPattern("*.gnu",dirName)
   copyPattern("*.m",dirName)


def runSimulation(newDir):
   currDir = os.getcwd()
   os.chdir(newDir)
   os.system("./setup initData")
   os.system("nohup nice -n 19 ./run initData 0 results &")
   os.chdir(currDir)

baseName = datetime.datetime.now().strftime("%y%m%d%H%M%S")+sys.argv[1];

re_s = [0.01]
res_s = [10,15,20,40,50]
porosity_s = [1.0]

for re in re_s:
   for res in res_s:
      for porosity in porosity_s:
         writeParametersMulti(porosity,res,re)
         compileProgram()
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/multi"+str(porosity)+"_"+str(res)+"_"+str(re)+"/"
         copyToDirectory(newDir)
         runSimulation(newDir)

