#include "customConstants.h"
#include "particle.h"
#include "io_data_vtk.h"
#include "misc.h"
#include "sph.h"
#include "sphIncompress.h"
#include "customOutput.h"
#include "dataLL.h"
#include <vector>
#include <ctime>
#include <iostream>
#include <fstream>

#define CHEB_N 65

int main(int argc, char *argv[]) {
   if (argc != 2) {
      cout << "Usage: setup outfilename" << endl;
      return(-1);
   }
   string filename = argv[1];
   
   CglobalVars globals;
   globals.mpiSize = 2;

   list<Cparticle> ps[globals.mpiSize];
   Cparticle p;
    
   p.tag = ps[0].size()+ps[1].size()+1;
   p.r = 0.2,0.2;
   p.dens = DENS;
   p.mass = PSEP*PSEP*DENS;
   p.h = H;
   p.v = 0,0;
   p.iam = sph;
   p.alpha = ALPHA;
   ps[0].push_back(p);

   p.tag = ps[0].size()+ps[1].size()+1;
   p.r = 0.9,0.9;
   p.dens = DENS;
   p.mass = PSEP*PSEP*DENS;
   p.h = H;
   p.v = 0,0;
   p.iam = sph;
   p.alpha = ALPHA;
   ps[1].push_back(p);

   p.tag = ps[0].size()+ps[1].size()+1;
   p.r = 0.2,0.4;
   p.dens = DENS;
   p.mass = PSEP*PSEP*DENS;
   p.h = H;
   p.v = VREF,0;
   p.iam = sph;
   p.alpha = ALPHA;
   ps[0].push_back(p);


   cout << "Total number of particles = " << ps[0].size()+ps[1].size() << endl;

 

   vector<vector<double> > vprocDomain(globals.mpiSize);
   vector<vector<int> > vprocNeighbrs(globals.mpiSize);

   vector<double> tmp(4);
   tmp[0] = 0;
   tmp[1] = 0.5;
   tmp[2] = 0;
   tmp[3] = 1;

   vprocDomain[0] = tmp;
   

   vprocDomain[1].resize(4);
   vprocDomain[1][0] = 0.5;
   vprocDomain[1][1] = 1;
   vprocDomain[1][2] = 0;
   vprocDomain[1][3] = 1;

   vprocNeighbrs[0].resize(4);
   vprocNeighbrs[0][0] = -1;
   vprocNeighbrs[0][1] = 1;
   vprocNeighbrs[0][2] = -1;
   vprocNeighbrs[0][3] = -1;

   vprocNeighbrs[1].resize(4);
   vprocNeighbrs[1][0] = 0;
   vprocNeighbrs[1][1] = -1;
   vprocNeighbrs[1][2] = -1;
   vprocNeighbrs[1][3] = -1;


   //cout << "creating dataLL"<<endl;
   //CdataLL *data = new CdataLL(ps,globals);
   //cout << "creating sph"<<endl;
   //CsphIncompress sph(data);
   //cout << "creating customOutput"<<endl;
   //CcustomOutput customOutput(data);
   //cout << "creating customSim"<<endl;
   //CcustomSim customSim(data,globals.time);


   cout << "Opening files for writing..."<<endl;
   Cio_data_vtk ioFile(filename.c_str());
   cout << "Writing Restart data to file..."<<endl;
   for (int i=0;i<globals.mpiSize;i++) {
      globals.mpiRank = i;
      ioFile.writeRestart(0,ps[i],&globals);
   }
   cout << "Writing Global data to file..."<<endl;
   ioFile.writeGlobals(0,&globals);
   cout << "Writing Domain data to file..."<<endl;
   ioFile.writeDomain(0,vprocDomain,vprocNeighbrs);
}

