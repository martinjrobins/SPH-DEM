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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define CHEB_N 65

int main(int argc, char *argv[]) {
   if (argc != 2) {
      cout << "Usage: setup outfilename" << endl;
      return(-1);
   }
   string filename = argv[1];

   vector<Cparticle> ps;
   Cparticle p;
   cout << "Creating box with sides: rmax = ["<<RMAX[0]<<" "<<RMAX[1]<<"] rmin = ["<<RMIN[0]<<" "<<RMIN[1]<<"]"<<endl;
   cout << "Reynolds Number = "<<REYNOLDS_NUMBER<<endl;
   cout << "Density = "<<DENS<<endl;
   cout << "number of particles on side = "<<NX<<endl;
   cout << "alpha = "<<ALPHA<<endl;
   cout << "viscosity = "<<VISCOSITY<<endl;
   cout << "maxtime = "<<MAXTIME<<endl;
   cout << "term vel (3d) is "<<(2.0/9.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,2)/(VISCOSITY*DENS)<<endl;
   cout << "term vel (2d) is "<<(1.0/3.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,2)/(VISCOSITY*DENS)<<endl;
 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<endl;

   p.tag = ps.size()+1;
   p.r[0] = (RMAX[0]-RMIN[0])/2;
   p.r[1] = (RMAX[1]-RMIN[1])/2;
   p.r[2] = (RMAX[2]-RMIN[2])/2;
   p.dens = DEM_DENS;
   p.h = DEM_RADIUS/2;
   p.mass = (4.0/3.0)*PI*pow(DEM_RADIUS,3)*DEM_DENS;
   p.v = 0.0;
   p.vhat = p.v;
   p.iam = dem;
   ps.push_back(p);

   const double totalDEMVol = p.mass/p.dens;
   //const double totalDEMVol = 0;
   const double newLiqDens = DENS*((RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0])*2*DEM_RADIUS-totalDEMVol)/((RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0])*2*DEM_RADIUS);
   cout << "after adding dem particles, new liquid density is "<<newLiqDens<<endl; 

   for (int i=0;i<NX;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=0;j<NY;j++) {
         for (int k=0;k<NZ;k++) {
            p.tag = ps.size()+1;
            p.r = (i+0.5)*PSEP+RMIN[0],(j+0.5)*PSEP+RMIN[1],(k+0.5)*PSEP+RMIN[2];
            p.dens = newLiqDens;
            p.mass = pow(PSEP,NDIM)*newLiqDens;
            p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
            p.v = 0.0;
            p.vhat = p.v;
            p.iam = sph;
            p.alpha = ALPHA;
            ps.push_back(p);
         }
      }
   }
   

   cout << "Total number of particles = " << ps.size() << endl;

   CglobalVars globals;

   vector<vector<double> > vprocDomain(globals.mpiSize);
   vector<Array<int,NDIM> > vprocNeighbrs(globals.mpiSize);
   vector<particleContainer > vps;
   vectInt split;
   split = 1,1,1;
   particleContainer pps;
   for (int i=0;i<ps.size();i++) {
      pps.push_back(ps[i]);
   }
   Nmisc::splitDomain(pps,split,vps,vprocDomain,vprocNeighbrs);


   cout << "Opening files for writing..."<<endl;
   Cio_data_vtk ioFile(filename.c_str(),&globals);
   cout << "Calculating Output stuff.."<<endl;
   //sph.calcOutputVars();
   //customOutput.calcOutput(0,&customSim,&ioFile);
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

