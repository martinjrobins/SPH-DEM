#!/usr/bin/python 
import sys
import os
import shutil,glob
import datetime
noFsolve = False
try:
   from scipy.optimize import fsolve
except ImportError:
   noFsolve = True
from numpy import exp,log10

def writePbsFile(name):
   filename = "pbs-script"
   file = open(filename,"w")
   file.write("#PBS -N "+name+"\n")
   file.write(open("pbs-script-base","r").read())
   file.close()

def writeParametersDry(res,dem_dens,dem_d):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double DEM_DENS = "+str(dem_dens)+";\n")
   file.write("const double DEM_RADIUS = "+str(dem_d/2.0)+";\n")
   file.close()

def writeParametersWet(res,dem_dens,dem_d):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define WET_START\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double DEM_DENS = "+str(dem_dens)+";\n")
   file.write("const double DEM_RADIUS = "+str(dem_d/2.0)+";\n")
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


def runSimulation(newDir,name):
   currDir = os.getcwd()
   os.chdir(newDir)
   os.system("./setup initData")
   shutil.copy(currDir+"/pbs-script-base",newDir)
   writePbsFile(name)
   os.system("qsub pbs-script")
   #os.system("nohup nice ./run initData 0 results &")
   os.chdir(currDir)



baseName = datetime.datetime.now().strftime("%y%m%d%H%M%S")+sys.argv[1];

dem_d = 1.1e-3;
dem_dens = 1160.0;

for res in [30]:
   writeParametersDry(res,dem_dens,dem_d)
   compileProgram()
   name = "dry"+str(res)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir)
   runSimulation(newDir,name)

for res in [30]:
   writeParametersWet(res,dem_dens,dem_d)
   compileProgram()
   name = "wet"+str(res)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir)
   runSimulation(newDir,name)

for res in [45,60]:
   writeParametersDry(res,dem_dens,dem_d)
   compileProgram()
   name = "dry"+str(res)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir)
   runSimulation(newDir,name)

for res in [45,60]:
   writeParametersWet(res,dem_dens,dem_d)
   compileProgram()
   name = "wet"+str(res)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir)
   runSimulation(newDir,name)


