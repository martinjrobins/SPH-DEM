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
 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<endl;

   ifstream fi("pos1.txt");
   if (fi.is_open()) {
      while (!fi.eof()) {
         string line;
         getline(fi,line);
         istringstream ss(line);
         ss >> p.r[0] >> p.r[1];
         p.tag = ps.size()+1;
         p.dens = DEM_DENS;
         p.h = DEM_RADIUS/2;
#ifdef _2D_DEM
         p.mass = PI*pow(DEM_RADIUS,2)*DEM_DENS;
#else
         p.mass = (4.0/3.0)*PI*pow(DEM_RADIUS,3)*DEM_DENS;
#endif
         p.v = 0.0;
         p.vhat = p.v;
         p.iam = dem;
         ps.push_back(p);
      }
   } else {
      cout << "Error reading pos1.txt" <<endl;
   }
   fi.close();
   

#ifdef _2D_DEM
   const double totalDEMVol = ps.size()*PI*pow(DEM_RADIUS,2);
   const double newLiqDens = DENS*((RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0])-totalDEMVol)/((RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0]));
   cout << "porosity = "<<1.0-totalDEMVol/((RMAX[0]-RMIN[0])*(RMAX[1]-RMIN[1]))<<endl;
#else
   const double totalDEMVol = ps.size()*(4.0/3.0)*PI*pow(DEM_RADIUS,3);
   const double newLiqDens = DENS*((RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0])*2*DEM_RADIUS-totalDEMVol)/((RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0])*2*DEM_RADIUS);
   cout << "porosity = "<<1.0-totalDEMVol/((RMAX[0]-RMIN[0])*(RMAX[1]-RMIN[1])*2*DEM_RADIUS)<<endl;
#endif
   //totalDEMVol = 0;
   cout << "after adding dem particles, new liquid density is "<<newLiqDens<<endl; 

   for (int i=0;i<NX;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=-2;j<=NY+2;j++) {
         if (j<=0) {
            p.iam = sphBoundary;
            p.dens = DENS;
            p.norm1[0] = 0;
            p.norm1[1] = j*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[1] = 1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if (j>=NY) {
            p.iam = sphBoundary;
            p.dens = DENS;
            p.norm1[0] = 0;
            p.norm1[1] = -j*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[1] = 1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else {
            p.iam = sph;
            p.dens = newLiqDens;
         }
         p.tag = ps.size()+1;
         p.r = (i+0.5)*PSEP+RMIN[0],j*PSEP+RMIN[1];
         p.mass = PSEP*PSEP*p.dens;
         p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
         p.v = 0.0;
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


   cout << "Opening files for writing..."<<endl;
   Cio_data_vtk ioFile(filename.c_str(),&globals);
   cout << "Calculating Output stuff.."<<endl;
   //sph.calcOutputVars();
   //customOutput.calcOutput(0,&customSim,&ioFile);
   cout << "Writing Restart data to file..."<<endl;
   int nProc = product(split);
   for (int i=0;i<nProc;i++) {
      globals.mpiRank = i;
      ioFile.setFilename(filename.c_str(),&globals);
      ioFile.writeGlobals(0,&globals);
      ioFile.writeRestart(0,vps[i],&globals);
      globals.mpiRank = 0;
   }
   cout << "Writing Global data to file..."<<endl;
   //ioFile.writeGlobals(0,&globals);
   ioFile.writeDomain(0,vprocDomain,vprocNeighbrs);

}

