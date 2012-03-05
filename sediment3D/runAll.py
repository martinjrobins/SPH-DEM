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

def writeParametersSingle(res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define LIQ_DEM_ONE_WAY_COUPLE\n")
   file.write("#define IN_WATER\n")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   #file.write("const int OUTSTEP = 200;");
   file.write("const double HMULT = "+str(2.0)+";\n")
   file.close()

def writeParametersSingleTwoWay(res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define IN_WATER\n")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(2.0)+";\n")
   file.close()

def writeParametersSingleTwoWayHigher(res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define IN_WATER\n")
   file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(60.0)+";\n")
   file.close()



def writeParametersSingleDiFelise(res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define IN_WATER\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(2.0)+";\n")
   file.close()

def writeParametersSingleDiFeliseHigher(res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define IN_WATER\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(60.0)+";\n")
   file.close()

def writeParametersSingleAM(res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define IN_WATER\n")
   #file.write("#define LIQ_DEM_ONE_WAY_COUPLE\n")
   #file.write("#define LIQ_DEM_SIMPLE_DRAG\n")
   file.write("#define LIQ_DEM_ADDED_MASS\n")
   file.write("const double POROSITY = 1.0;\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(2.0)+";\n")
   file.close()

def writeParametersMulti(p,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define IN_WATER\n")
   file.write("#define MANY_PARTICLES\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const double DEM_DENS = DENS*2.5;\n")
   file.write("const int NX = 10;\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(2.0)+";\n")
   file.close()


def writeParametersMultiGrid(p,res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define MANY_PARTICLES\n")
   file.write("#define GRID_OF_DEM\n")
   file.write("#define GRID_OF_DEM_PERT\n")
   file.write("#define IN_WATER\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(2.0)+";\n")
   file.close()

def writeParametersMultiGridShort(p,res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define MANY_PARTICLES\n")
   file.write("#define GRID_OF_DEM\n")
   file.write("#define GRID_OF_DEM_PERT\n")
   file.write("#define IN_WATER\n")
   file.write("#define SHORT_TIME\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(2.0)+";\n")
   file.close()



def writeParametersMultiGridBlock(p,res,dens,visc,re):
   filename = "parameters.h"
   file = open(filename,"w")
   file.write("#define IN_WATER\n")
   file.write("#define MANY_PARTICLES\n")
   file.write("#define GRID_OF_DEM\n")
   file.write("#define DEM_BLOCK\n")
   file.write("const double POROSITY = "+str(p)+";\n")
   file.write("const int NX = "+str(res)+";\n")
   file.write("const double VISCOSITY = "+str(visc)+";\n")
   file.write("const double DENS = "+str(dens)+";\n")
   file.write("const double REYNOLDS_NUMBER = "+str(re)+";\n")
   file.write("const double HMULT  = "+str(2.0)+";\n")
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
   #shutil.copy(currDir+"/pbs-script-base",newDir)
   #writePbsFile(name)
   #os.system("qsub pbs-script")
   #os.system("nohup nice -n 19 ./run initData 0 results &")
   os.chdir(currDir)

class fluid:
   pass

def f(re,args):
   beta = 3.7 - 0.65*exp(-((1.5-log10(re))**2.0)/2.0);
   factor = args[1]**(-beta);
   f = 0.3969*re**2 + 6.048*re**1.5 + 23.04*re - (4.0/3.0)*args[0]/factor;
   return f

def calcRe(d,dem_dens,dens,visc,p):
   ar = d**3*(dem_dens/dens - 1.0)*9.81/(visc**2);
   re0 = ar/18.0;
   global noFsolve
   if not noFsolve:
      rePlus = fsolve(f,re0,[ar,p]);
   else:
      rePlus = re0;
   re = rePlus;
   return re

def f2(re,args):
   f = re + 0.2652*re**1.5 + 0.172*re**2 - args[0]/18.0;
   return f

def calcRe2(d,dem_dens,dens,visc):
   ar = d**3*(dem_dens/dens - 1.0)*9.81/(visc**2);
   re0 = ar/18.0;
   global noFsolve
   if not noFsolve:
      rePlus = fsolve(f,re0,[ar,p]);
   else:
      rePlus = re0;
   re = rePlus;
   return re



baseName = datetime.datetime.now().strftime("%y%m%d%H%M%S")+sys.argv[1];

water_glycerol = fluid();
water_glycerol.dens = 1150.0;
water_glycerol.visc = 8.9e-3/water_glycerol.dens;
water_glycerol.name = "water_glycerol";
water = fluid();
water.dens = 1000.0;
water.visc = 8.9e-4/water.dens;
water.name = "water";
water2 = fluid();
water2.dens = 1000.0;
water2.visc = 0.6*8.9e-4/water.dens;
water2.name = "water2Re2Re";
water3 = fluid();
water3.dens = 1000.0;
water3.visc = 0.24*8.9e-4/water.dens;
water3.name = "water10Re";
air = fluid();
air.dens = 1.1839;
air.visc = 18.6e-6/air.dens;
air.name = "air";

d = 100.0e-6;
dem_dens = 2500.0;

for i in [water_glycerol,water,air,water2,water3]:
   #re = calcRe2(d,dem_dens,i.dens,i.visc);
   #termV = re*i.visc/(d);
   #print "dens = ",i.dens," kg/m^3 kinematic viscosity = ",i.visc," m^2/s porosity = ",1.0," RE_p = ",re," term velocity = ",termV," m/s";
   for p in [0.5,0.55,0.6,0.65,0.7,0.8,0.86,0.9,1.0]:
      re = calcRe(d,dem_dens,i.dens,i.visc,p);

      beta = 3.7 - 0.65*exp(-((1.5-log10(re))**2.0)/2.0);
      C =(0.63 + 4.8/(re**0.5))**2.0;
      termV = re*i.visc/(d);
      print "dens = ",i.dens," kg/m^3 kinematic viscosity = ",i.visc," m^2/s porosity = ",p," RE_p = ",re," term velocity = ",termV," m/s"," beta = ",beta," f = ",C*p**(2-beta)," Ar = ",d**3*i.dens*(dem_dens-i.dens)*9.81/((i.visc*i.dens)**2);

for i in [water_glycerol,water,air]:
   re = calcRe(d,dem_dens,i.dens,i.visc,1.0);
   writeParametersSingle(10,i.dens,i.visc,re)
   compileProgram()
   name = "single"+str(10)+"_"+i.name+"/"
   newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
   copyToDirectory(newDir)
   runSimulation(newDir,name)


for i in [water]:
   for res in [10,20,30,40]:
      re = calcRe(d,dem_dens,i.dens,i.visc,1.0);
      writeParametersSingleTwoWay(res,i.dens,i.visc,re)
      compileProgram()
      name = "singleTwoWay"+str(res)+"_"+i.name+"/"
      newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
      copyToDirectory(newDir)
      runSimulation(newDir,name)

for i in [air,water_glycerol,water2,water3]:
   for res in [10]:
      re = calcRe(d,dem_dens,i.dens,i.visc,1.0);
      writeParametersSingleTwoWay(res,i.dens,i.visc,re)
      compileProgram()
      name = "singleTwoWay"+str(res)+"_"+i.name+"/"
      newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
      copyToDirectory(newDir)
      runSimulation(newDir,name)

for i in [air]:
   for res in [10]:
      re = calcRe(d,dem_dens,i.dens,i.visc,1.0);
      writeParametersSingleTwoWayHigher(res,i.dens,i.visc,re)
      compileProgram()
      name = "singleTwoWay"+str(res)+"Higher_"+i.name+"/"
      newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
      copyToDirectory(newDir)
      runSimulation(newDir,name)
      
      re = calcRe(d,dem_dens,i.dens,i.visc,1.0);
      writeParametersSingleDiFeliseHigher(res,i.dens,i.visc,re)
      compileProgram()
      name = "singleDiFelise"+str(res)+"Higher_"+i.name+"/"
      newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
      copyToDirectory(newDir)
      runSimulation(newDir,name)

for i in [water_glycerol,water,air,water2,water3]:
   for res in [10]:
      re = calcRe(d,dem_dens,i.dens,i.visc,1.0);
      writeParametersSingleDiFelise(res,i.dens,i.visc,re)
      compileProgram()
      name = "singleDiFelise"+str(res)+"_"+i.name+"/"
      newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
      copyToDirectory(newDir)
      runSimulation(newDir,name)

for i in [water_glycerol,water,air]:
   for res in [10]:
      re = calcRe(d,dem_dens,i.dens,i.visc,1.0);
      writeParametersSingleAM(res,i.dens,i.visc,re)
      compileProgram()
      name = "singleAM"+str(res)+"_"+i.name+"/"
      newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
      copyToDirectory(newDir)
      runSimulation(newDir,name)

for i in [water]:
   for res in [10,20,30,40]:
      for porosity in [0.8]:
         re = calcRe(d,dem_dens,i.dens,i.visc,porosity);
         writeParametersMultiGridBlock(porosity,res,i.dens,i.visc,re)
         compileProgram()
         name = "multi"+str(porosity)+"gridBlock"+str(res)+"_"+i.name+"/"
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
         copyToDirectory(newDir)
         runSimulation(newDir,name)

for i in [water_glycerol,water,air]:
   for res in [10]:
      for porosity in [0.6,0.7,0.8,0.9,0.95]:
         re = calcRe(d,dem_dens,i.dens,i.visc,porosity);
         writeParametersMultiGridBlock(porosity,res,i.dens,i.visc,re)
         compileProgram()
         name = "multi"+str(porosity)+"gridBlock"+str(res)+"_"+i.name+"/"
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
         copyToDirectory(newDir)
         runSimulation(newDir,name)

for i in [water_glycerol,water,air]:
   for res in [10]:
      for porosity in [0.8]:
         re = calcRe(d,dem_dens,i.dens,i.visc,porosity);
         writeParametersMultiGrid(porosity,res,i.dens,i.visc,re)
         compileProgram()
         name = "multi"+str(porosity)+"grid"+str(res)+"_"+i.name+"/"
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
         copyToDirectory(newDir)
         runSimulation(newDir,name)

for i in [water]:
   for res in [12,15]:
      for porosity in [0.8]:
         re = calcRe(d,dem_dens,i.dens,i.visc,porosity);
         writeParametersMultiGridBlock(porosity,res,i.dens,i.visc,re)
         compileProgram()
         name = "multi"+str(porosity)+"gridBlock"+str(res)+"_"+i.name+"/"
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
         copyToDirectory(newDir)
         runSimulation(newDir,name)

for i in [water_glycerol,water]:
   for res in [10]:
      for porosity in [0.5,0.55,0.65]:
         re = calcRe(d,dem_dens,i.dens,i.visc,porosity);
         writeParametersMultiGridBlock(porosity,res,i.dens,i.visc,re)
         compileProgram()
         name = "multi"+str(porosity)+"gridBlock"+str(res)+"_"+i.name+"/"
         newDir = os.environ["HOME"]+"/data/sediment3D/"+baseName+"/"+name
         copyToDirectory(newDir)
         runSimulation(newDir,name)


