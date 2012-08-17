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

def writePbsFile(name,ncpu,restart):
   filename = "pbs-script"
   file = open(filename,"w")
   file.write("#PBS -N "+name+"\n")  
   file.write("#PBS -l nodes=1:ppn="+str(ncpu)+"\n")
   file.write(open("pbs-script-base","r").read())
   if restart:
      file.write("mpiexec -np "+str(ncpu)+" ./run results 99 results\n")
   else:
      file.write("mpiexec -np "+str(ncpu)+" ./run initData 0 results\n")
   file.close()

def writeParametersDry(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define NIX "+str(res)+"\n")
   file.write("const double DEM_DENS = "+str(dem_dens)+";\n")
   file.write("const double DEM_RADIUS = "+str(dem_d/2.0)+";\n")
   #file.write("const double INLET_RADIUS = "+str(ri)+";\n")
   file.write("const double REAL_INLET_RADIUS = "+str(rri)+";\n")
   file.write("const double INLET_FLOW_RATE = "+str(vi)+";\n")
   file.write("#define LINEAR\n")
   file.write("const double NCPU_X = "+str(ncpu_x)+";\n")
   file.write("const double NCPU_Y = "+str(ncpu_y)+";\n")
   file.write("const double NCPU_Z = "+str(ncpu_z)+";\n")
   file.close()

def writeParametersDryLub(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define NIX "+str(res)+"\n")
   file.write("const double DEM_DENS = "+str(dem_dens)+";\n")
   file.write("const double DEM_RADIUS = "+str(dem_d/2.0)+";\n")
   #file.write("const double INLET_RADIUS = "+str(ri)+";\n")
   file.write("const double REAL_INLET_RADIUS = "+str(rri)+";\n")
   file.write("const double INLET_FLOW_RATE = "+str(vi)+";\n")
   file.write("#define LUBRICATION\n")
   file.write("const double NCPU_X = "+str(ncpu_x)+";\n")
   file.write("const double NCPU_Y = "+str(ncpu_y)+";\n")
   file.write("const double NCPU_Z = "+str(ncpu_z)+";\n")
   file.close()


def writeParametersWet(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define NIX "+str(res)+"\n")
   file.write("#define WET_START\n")
   file.write("const double DEM_DENS = "+str(dem_dens)+";\n")
   file.write("const double DEM_RADIUS = "+str(dem_d/2.0)+";\n")
   #file.write("const double INLET_RADIUS = "+str(ri)+";\n")
   file.write("const double REAL_INLET_RADIUS = "+str(rri)+";\n")
   file.write("const double INLET_FLOW_RATE = "+str(vi)+";\n")
   file.write("#define LINEAR\n")
   file.write("const double NCPU_X = "+str(ncpu_x)+";\n")
   file.write("const double NCPU_Y = "+str(ncpu_y)+";\n")
   file.write("const double NCPU_Z = "+str(ncpu_z)+";\n")
   file.close()

def writeParametersWetLub(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define NIX "+str(res)+"\n")
   file.write("#define WET_START\n")
   file.write("const double DEM_DENS = "+str(dem_dens)+";\n")
   file.write("const double DEM_RADIUS = "+str(dem_d/2.0)+";\n")
   #file.write("const double INLET_RADIUS = "+str(ri)+";\n")
   file.write("const double REAL_INLET_RADIUS = "+str(rri)+";\n")
   file.write("const double INLET_FLOW_RATE = "+str(vi)+";\n")
   file.write("#define LUBRICATION\n")
   file.write("const double NCPU_X = "+str(ncpu_x)+";\n")
   file.write("const double NCPU_Y = "+str(ncpu_y)+";\n")
   file.write("const double NCPU_Z = "+str(ncpu_z)+";\n")
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

def copyToDirectory(dirName,dirRestart):
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
   copyPattern(dirRestart+"/resultsRestart0000099*.*",dirName)
   copyPattern(dirRestart+"/resultsGlobals*",dirName)
   copyPattern(dirRestart+"/resultsDomain*",dirName)


def runSimulation(newDir,name,ncpu,restart):
   currDir = os.getcwd()
   os.chdir(newDir) 
   shutil.copy(currDir+"/pbs-script-base",newDir)
   writePbsFile(name,ncpu,restart)
   if not restart:
      os.system("./setup initData")
   os.system("qsub pbs-script")
   os.chdir(currDir)



baseName = datetime.datetime.now().strftime("%y%m%d%H%M%S")+sys.argv[1];

dem_d = 1.1e-3;
dem_dens = 1160.0;


ncpu_x = 1;
ncpu_y = 1;
ncpu_z = 1;
ncpu = ncpu_x*ncpu_y*ncpu_z;

#for res in [60]:
#   writeParametersDry(res,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
#   compileProgram()
#   name = "dry"+str(res)
#   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
#   copyToDirectory(newDir)
#   runSimulation(newDir,name,ncpu)


#for res in [60]:
#   writeParametersDryLub(res,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
#   compileProgram()
#   name = "dryLubrication"+str(res)
#   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
#   copyToDirectory(newDir)
#   runSimulation(newDir,name,ncpu)

#for res in [60]:
#   writeParametersWet(res,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
#   compileProgram()
#   name = "wet"+str(res)
#   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
#   copyToDirectory(newDir)
#   runSimulation(newDir,name,ncpu)

ncpu_x = 1;
ncpu_y = 1;
ncpu_z = 1;
ncpu = ncpu_x*ncpu_y*ncpu_z;

ri = 1.0/1000.0;
#ri = 3.275/1000.0;
#rri = 2.5/1000.0;
rri = 1.0/1000.0;
vi = 100.0

#dirRestart = "/home/mrobinson/data/cell3D/120606145231varResBaseParam"
dirRestart = "."
restart = False



ncpu_x = 2;
ncpu_y = 1;
ncpu_z = 1;
ncpu = ncpu_x*ncpu_y*ncpu_z;

for res in [1]:
   writeParametersWet(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
   compileProgram()
   name = "wet"+str(res)+"cpu"+str(ncpu)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir,dirRestart)
   runSimulation(newDir,name,ncpu,restart)

   writeParametersDry(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
   compileProgram()
   name = "dry"+str(res)+"cpu"+str(ncpu)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir,dirRestart)
   runSimulation(newDir,name,ncpu,restart)



ncpu_x = 2;
ncpu_y = 2;
ncpu_z = 2;
ncpu = ncpu_x*ncpu_y*ncpu_z;

for res in [2]:
   writeParametersWet(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
   compileProgram()
   name = "wet"+str(res)+"cpu"+str(ncpu)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir,dirRestart)
   runSimulation(newDir,name,ncpu,restart)

   writeParametersDry(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
   compileProgram()
   name = "dry"+str(res)+"cpu"+str(ncpu)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir,dirRestart)
   runSimulation(newDir,name,ncpu,restart)

ncpu_x = 2;
ncpu_y = 2;
ncpu_z = 3;
ncpu = ncpu_x*ncpu_y*ncpu_z;

for res in [3,4]:
   writeParametersWet(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
   compileProgram()
   name = "wet"+str(res)+"cpu"+str(ncpu)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir,dirRestart)
   runSimulation(newDir,name,ncpu,restart)

   writeParametersDry(res,ri,rri,vi,dem_dens,dem_d,ncpu_x,ncpu_y,ncpu_z)
   compileProgram()
   name = "dry"+str(res)+"cpu"+str(ncpu)
   newDir = os.environ["HOME"]+"/data/cell3D/"+baseName+"/"+name+"/"
   copyToDirectory(newDir,dirRestart)
   runSimulation(newDir,name,ncpu,restart)

