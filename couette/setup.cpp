#include "customConstants.h"
#include "io_data_vtk.h"
#include "misc.h"
#include "sphIncompress.h"
#include "customOutput.h"
#include "dataLL.h"
#include "particle.h"
#include "sph.h"
#include <vector>
#include <iostream>

int main(int argc, char *argv[]) {
   if (argc != 2) {
      cout << "Usage: setup outfilename" << endl;
      return(-1);
   }
   string filename = argv[1];



   particleContainer ps;
   Cparticle p;
   cout << "Creating box with sides: rmax = ["<<RMAX[0]<<" "<<RMAX[1]<<"] rmin = ["<<RMIN[0]<<" "<<RMIN[1]<<"]"<<endl;
   cout << "Reynolds Number = "<<REYNOLDS_NUMBER<<endl;
   cout << "Density = "<<DENS<<endl;
   cout << "number of particles on side = "<<NX<<endl;
   cout << "alpha = "<<ALPHA<<endl;
   cout << "viscosity = "<<VISCOSITY<<endl;
   cout << "maxtime = "<<MAXTIME<<endl;


   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<endl;
   for (int i=0;i<NX;i++) {
      cout<<"i = "<<i<<endl;
      for (int j=0;j<=NY;j++) {
         p.r = (i+0.5)*PSEP+RMIN[0],j*PSEP+RMIN[1];
         p.dens = DENS;
         p.mass = PSEP*PSEP*DENS;
         p.h = H;
         if (j<3) {
            p.v = 0,0;
            p.iam = sphBoundary;
         } else if (j>NY-3) {
            p.v = VREF,0;
            p.v[0] += (j-(NY-2))*VREF/(NY-5);
            p.iam = sphBoundary;
         } else {
            p.v = 0,0;
            p.iam = sph;
         }
         p.vhat = p.v;
         p.alpha = ALPHA;
         ps.push_back(p);
      }
   }

   cout << "Total number of particles = " << ps.size() << endl;
   CglobalVars globals;
 
   vector<vector<double> > vprocDomain(globals.mpiSize);
   vector<Array<int,NDIM> > vprocNeighbrs(globals.mpiSize);
   vector<particleContainer > vps;
   vectInt split;
   split = 1,1;
   
   particleContainer pps;
   for (int i=0;i<ps.size();i++) {
      pps.push_back(ps[i]);
   }
   Nmisc::splitDomain(pps,split,vps,vprocDomain,vprocNeighbrs);


   
   //CdataLL *data = new CdataLL(ps,globals);
   //CsphIncompress sph(data);
   //CcustomOutput customOutput(data);
   //CcustomSim customSim(data,globals.time);

   //cout << "Calculating Output stuff.."<<endl;
   //sph.calcOutputVars();

   cout << "Opening files for writing..."<<endl;
   Cio_data_vtk ioFile(filename.c_str(),&globals);
   cout << "Writing Restart data to file..."<<endl;
   int nProc = product(split);
   for (int i=0;i<nProc;i++) {
      globals.mpiRank = i;
      for (int j=0;j<NDIM*2;j++)
         globals.procDomain[j] = vprocDomain[i][j];
      globals.procNeighbrs = vprocNeighbrs[i];
      ioFile.setFilename(filename.c_str(),&globals);
      ioFile.writeGlobals(0,&globals);
      ioFile.writeRestart(0,vps[i],&globals);
      ioFile.writeDomain(0,&globals);
      globals.mpiRank = 0;
   }

}

