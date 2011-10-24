#! /usr/bin/env python
import glob
import string
import os
import math


listFiles = glob.glob('*')
dirs = []
for filename in listFiles:
   if os.path.isdir(filename):
      dirs.append(filename)

n = len(dirs)
numCPUs = 2
nperCPU = int(math.ceil(n/numCPUs))
timePerRun = 0.2
totalTime = int(math.ceil(0.2*nperCPU))

print "Running ",n," simulations over ",numCPUs," CPUS. Total time = ",totalTime

pbsScripts = []
for i in range(0,numCPUs):
   file = open("pbsScript"+str(i),"w")
   file.write("#!/bin/bash\n")
   file.write("#PBS -N runParams"+str(i)+"\n")
   file.write("#PBS -l nodes=1\n")
   file.write("#PBS -l walltime="+str(totalTime*2)+":00:00\n")
   file.write("cd $PBS_O_WORKDIR\n")
   pbsScripts.append(file)


i = 0
for dir in dirs:
   cpu = i%numCPUs
   (w,v,h) = string.split(dir,"-")

   print "Processing (w,v,h) = ",w,v,h

   file = pbsScripts[cpu]
   file.write("cd "+dir+"\n")
   file.write("./setup initData\n")
   file.write("./run initData 0 results\n")
   #file.write("./calcLyap results 0 2.33123\n")
   file.write("cd ..\n")
   i = i + 1

for i in range(0,numCPUs):
   file = pbsScripts[i]
   file.close()

