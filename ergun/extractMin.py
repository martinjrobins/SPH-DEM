import csv
from string import atof

dataCol = 19
data = csv.reader(open('demVel.dat','rb'),delimiter=' ')
dataWrite = csv.writer(open('minDem.dat', 'wb'), delimiter=' ')
first = data.next()
first = data.next()
time = atof(first[0])
minn = atof(first[dataCol])
for row in data:
   cTime = atof(row[0])
   cData = atof(row[dataCol])
   if (time != cTime):
      dataWrite.writerow([time,minn])
      time = cTime
   else:
      if (minn > cData):
         minn = cData


