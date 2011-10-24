import vtk
import glob
import os
import numpy
import sys

def addArray(ug,array,name):
   theArray = vtk.vtkFloatArray()
   theArray.SetName(name)
   #theArray.SetNumberOfValues(len(array))
   #for i in array:
   #   theArray.InsertNextValue(i)
   theArray.SetArray(array,len(array),1)
   ug.GetPointData.AddArray(theArray)

def getStr(data,T):
   n = data.GetNumberOfTuples()
   str = numpy.zeros(n)
   for i in range(0,n):
      str[i] = exp(T*data.GetValue(i))
   return str

def getFTLE(data):
   n = data.GetNumberOfTuples()
   str = numpy.zeros(n)
   for i in range(0,n):
      str[i] = data.GetValue(i)
   return str

try:
   theFilename = sys.argv[1]
   T = atof(sys.argv[1])
except IndexError:
   print 'usage: %s filename T'%sys.argv[0]

listFiles = glob.glob('*')
dirs = []
for filename in listFiles:
   if os.path.isdir(filename):
      dirs.append(filename)

reader = vtk.vtkUnstructuredGridReader()
ugRead = vtk.vtkUnstructuredGrid()

FaveStr = []
FstdDevStr = []
FmedianStr = []
BaveStr = []
BstdDevStr = []
BmedianStr = []
FaveFTLE = []
FstdDevFTLE = []
FmedianFTLE = []
BaveFTLE = []
BstdDevFTLE = []
BmedianFTLE = []
for dir in dirs:
   reader.setFileName(dir+"/"+theFilename)
   reader.SetOutput(ugRead)
   reader.Update()
   n = ugRead.GetNumberOfPoints()
   data = ugRead.GetPointData().GetArray("ftle")
   dataStr = getStr(data,T)
   dataFTLE = getFTLE(data)

   FaveStr.append(numpy.mean(dataStr))
   FstdDevStr.append(numpy.std(dataStr))
   FmedianStr.append(numpy.median(dataStr))
   FaveFTLE.append(numpy.mean(dataFTLE))
   FstdDevFTLE.append(numpy.std(dataFTLE))
   FmedianFTLE.append(numpy.median(dataFTLE))

   reader.setFileName(dir+"/"+ftleBfilename)
   reader.SetOutput(ugRead)
   reader.Update()
   n = ugRead.GetNumberOfPoints()
   data = ugRead.GetPointData().GetArray("ftle")

   dataStr = getStr(data,T)
   dataFTLE = getFTLE(data)

   BaveStr.append(numpy.mean(dataStr))
   BstdDevStr.append(numpy.std(dataStr))
   BmedianStr.append(numpy.median(dataStr))
   BaveFTLE.append(numpy.mean(dataFTLE))
   BstdDevFTLE.append(numpy.std(dataFTLE))
   BmedianFTLE.append(numpy.median(dataFTLE))

ugNew = vtk.vtkUnstructuredGrid()
ugNew.SetPoints(ugRead.GetPoints())

addArray(points,FaveStr,"FaveStr")
addArray(points,FstdDevStr,"FstdDevStr")
addArray(points,FmedianStr,"FmedianStr")
addArray(points,FaveFTLE,"FaveFTLE")
addArray(points,FstdDevFTLE,"FstdDevFTLE")
addArray(points,FmedianFTLE,"FmedianFTLE")
addArray(points,BaveStr,"BaveStr")
addArray(points,BstdDevStr,"BstdDevStr")
addArray(points,BmedianStr,"BmedianStr")
addArray(points,BaveFTLE,"BaveFTLE")
addArray(points,BstdDevFTLE,"BstdDevFTLE")
addArray(points,BmedianFTLE,"BmedianFTLE")

writer = vtkXMLUnstructuredGridWriter()
writer.SetInput(ugNew)
writer.SetDataModeToBinary()
writer.SetFileName("paramsResult.vtu")
writer.Write()

