#include <fstream>
#include <iostream>
#include "customConstants.h"
#include <string>
#include <iostream>
#include <sstream>
using namespace std;

int main(int argc, char *argv[]) {
   string filename = argv[1];
   string filenameIn = filename + "Globals0.dat";
   string filenameOut = filename + ".pvd";

   ofstream fo(filenameOut.c_str());
   ifstream fi(filenameIn.c_str());
   
   string line;
   getline(fi,line);

   fo << "<?xml version=\"1.0\"?>"<<endl;
   fo << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" <<endl;
   fo << "<Collection>" <<endl;

   int currstep,sphStep,oldStep;
   double time;
   currstep = 0;
   for (int i=0;i<OUTSTEP;i++) {
      getline(fi,line);
      istringstream ss(line);
      oldStep = currstep;
      ss >> currstep >> sphStep >> time;
      if (fi.eof()||(currstep != oldStep+1)) break;

      char strTimestep[20];
      sprintf(strTimestep,"%7.7d",currstep);
      string filenameStep = filename+strTimestep+".pvtu";

      fo << "<DataSet timestep=\"" << time << "\" group=\"\" part = \"0\" file=\"" << filenameStep << "\"/>" <<endl;
   }
   fo << "</Collection>" <<endl;
   fo << "</VTKFile>" <<endl;

   fo.close();
   fi.close();
}
