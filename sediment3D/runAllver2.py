#!/usr/bin/python 
import sys
import os
import shutil,glob
import datetime

def writeParametersSingle(res,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define LIQ_DEM_ONE_WAY_COUPLE\n")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   #file.write("const int OUTSTEP = 200;");
   file.close()

def writeParametersSingleTwoWay(res,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.close()

def writeParametersSingleDiFelise(res,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.close()


def writeParametersSingleAM(re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define LIQ_DEM_ONE_WAY_COUPLE\n")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("#define LIQ_DEM_ADDED_MASS\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = 10;\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.close()

def writeParametersMulti(p):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define MANY_PARTICLES\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const double DEM_DENS = DENS*8.0;\n")
   file.write("const int NX = 10;\n")
   file.write("const double REYNOLDS_NUMBER = 0.01;\n")
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

def writeParametersMultiGridBlock(p,res,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define MANY_PARTICLES\n")
   file.write("#define GRID_OF_DEM\n")
   file.write("#define DEM_BLOCK\n")
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
   #os.system("nohup nice -n 19 ./run initData 0 results &")
   os.chdir(currDir)

baseName = datetime.datetime.now().strftime("%y%m%d%H%M%S")+sys.argv[1];

writeParametersSingle(10,0.01)
compileProgram()
newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/single"+str(10)+"_"+str(0.01)+"/"
copyToDirectory(newDir)
runSimulation(newDir)

re_s = [0.01]
res_s = [10,20,30,40,50,60]
porosity_s = [1.0]

for re in re_s:
   for res in res_s:
      for porosity in porosity_s:
         writeParametersSingleTwoWay(res,re)
         compileProgram()
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/singleTwoWay"+str(res)+"_"+str(re)+"/"
         copyToDirectory(newDir)
         runSimulation(newDir)

re_s = [0.1,1,20,50,100]
res_s = [10]
porosity_s = [1.0]

for re in re_s:
   for res in res_s:
      for porosity in porosity_s:
         writeParametersSingleDiFelise(res,re)
         compileProgram()
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/singleDiFelise"+str(res)+"_"+str(re)+"/"
         copyToDirectory(newDir)
         runSimulation(newDir)


re_s = [0.01]
res_s = [10,20,30,40,50,60]
porosity_s = [0.8]

for re in re_s:
   for res in res_s:
      for porosity in porosity_s:
         writeParametersMultiGridBlock(porosity,res,re)
         compileProgram()
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/multi"+str(porosity)+"gridBlock"+str(res)+"_"+str(re)+"/"
         copyToDirectory(newDir)
         runSimulation(newDir)

re_s = [0.01,1.0,50.0]
res_s = [10]
porosity_s = [0.45,0.5,0.55,0.6,0.7,0.8,0.9,0.95,0.99]

for re in re_s:
   for res in res_s:
      for porosity in porosity_s:
         writeParametersMultiGridBlock(porosity,res,re)
         compileProgram()
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/multi"+str(porosity)+"gridBlock"+str(res)+"_"+str(re)+"/"
         copyToDirectory(newDir)
         runSimulation(newDir)

re_s = [0.01,1.0]
res_s = [10]
porosity_s = [0.8]

for re in re_s:
   for res in res_s:
      for porosity in porosity_s:
         writeParametersMultiGrid(porosity,res,re)
         compileProgram()
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/multi"+str(porosity)+"grid"+str(res)+"_"+str(re)+"/"
         copyToDirectory(newDir)
         runSimulation(newDir)
