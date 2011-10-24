#!/usr/bin/python 
import sys
import os
import shutil
import datetime

def writeParametersSingle(res):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define LIQ_DEM_ONE_WAY_COUPLE\n")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const double DEM_DENS = DENS*8.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.close()

def writeParametersSingleTwoWay(res):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const double DEM_DENS = DENS*8.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.close()


def writeParametersSingleAM():
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define LIQ_DEM_ONE_WAY_COUPLE\n")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("#define LIQ_DEM_ADDED_MASS\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const double DEM_DENS = DENS*8.0;\n")
   file.write("const int NX = 10;\n")
   file.close()

def writeParametersMulti(p):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define MANY_PARTICLES\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const double DEM_DENS = DENS*2.0;\n")
   file.write("const int NX = 10;\n")
   file.close()


def writeParametersMultiGrid(p):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define MANY_PARTICLES\n")
   file.write("#define GRID_OF_DEM\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const double DEM_DENS = DENS*2.0;\n")
   file.write("const int NX = 10;\n")
   file.close()



def myRemove(name):
   if os.path.exists(name):
      os.remove(name)

def compileProgram():
   os.system("make clean")
   os.system("make -j 8")

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


def runSimulation(newDir):
   currDir = os.getcwd()
   os.chdir(newDir)
   os.system("./setup initData")
   os.system("nohup nice ./run initData 0 results &")
   os.chdir(currDir)

baseName = datetime.datetime.now().strftime("%y%m%d%H%M%S")+sys.argv[1];

#writeParametersSingle(10)
#compileProgram()
#newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/single/"
#copyToDirectory(newDir)
#runSimulation(newDir)

#writeParametersSingleTwoWay(10)
#compileProgram()
#newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/single10/"
#copyToDirectory(newDir)
#runSimulation(newDir)

#writeParametersSingleTwoWay(15)
#compileProgram()
#newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/single15/"
#copyToDirectory(newDir)
#runSimulation(newDir)

#writeParametersSingleTwoWay(20)
#compileProgram()
#newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/single20/"
#copyToDirectory(newDir)
#runSimulation(newDir)


#writeParametersSingleAM()
#compileProgram()
#newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/singleAM/"
#copyToDirectory(newDir)
#runSimulation(newDir)

#writeParametersMulti(0.7)
#compileProgram()
#newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/multi70/"
#copyToDirectory(newDir)
#runSimulation(newDir)

#writeParametersMulti(0.8)
#compileProgram()
#newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/multi80/"
#copyToDirectory(newDir)
#runSimulation(newDir)

#writeParametersMulti(0.9)
#compileProgram()
#newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/multi90/"
#copyToDirectory(newDir)
#runSimulation(newDir)

writeParametersMultiGrid(0.9)
compileProgram()
newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/multi90grid/"
copyToDirectory(newDir)
runSimulation(newDir)

