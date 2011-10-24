import os
import shutil

def writeParameters(w,v,h):
   filename = "mixerParameters.h"
   file = open(filename,"w")
   file.write("const double WRATIO = "+str(w)+";\n")
   file.write("const double VRATIO = "+str(v)+";\n")
   file.write("const double HRATIO = "+str(h)+";\n")
   file.close()

def myRemove(name):
   if os.path.exists(name):
      os.remove(name)

def compileProgram():
   myRemove("setup.o")
   myRemove("customSim.o")
   myRemove("customOutput.o")
   os.system("make")
   os.system("make post")
   os.system("make calcLyap")

def copyToDirectory(w,v,h,baseDir):
   dirName = baseDir+"/"+str(w)+"-"+str(v)+"-"+str(h)
   d = os.path.dirname(dirName);
   if not os.path.exists(dirName):
      os.makedirs(dirName)
   dirName = dirName+"/"
   shutil.copy("run",dirName)
   shutil.copy("setup",dirName)
   shutil.copy("calcLyap",dirName)
   shutil.copy("post",dirName)
   shutil.copy("createPVD",dirName)
   shutil.copy("customConstants.h",dirName)
   shutil.copy("mixerParameters.h",dirName)


wratios = [2.0/5.0, 1.0]
vratios = [0.0,1.0]
hratios = [0.0,2.0]

for w in wratios:
   for v in vratios:
      for h in hratios:
         print "\ngenerating (w,v,h) = ",w,v,h,"\n"
         writeParameters(w,v,h)
         compileProgram()
         copyToDirectory(w,v,h,os.environ["HOME"]+"/data/2dMicro/paramSearch")

