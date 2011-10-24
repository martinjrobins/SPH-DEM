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
   cout << "Creating box with sides: rmax = ["<<BMAX[0]<<" "<<BMAX[1]<<"] rmin = ["<<BMIN[0]<<" "<<BMIN[1]<<"]"<<endl;
   cout << "Reynolds Number = "<<REYNOLDS_NUMBER<<endl;
   cout << "Density = "<<DENS<<endl;
   cout << "number of particles on side = "<<NX<<endl;
   cout << "alpha = "<<ALPHA<<endl;
   cout << "viscosity = "<<VISCOSITY<<endl;
   cout << "maxtime = "<<MAXTIME<<endl;
   cout << "term vel (3d) is "<<(2.0/9.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,2)/(VISCOSITY*DENS)<<endl;
   cout << "term vel (2d) is "<<(1.0/3.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,2)/(VISCOSITY*DENS)<<endl;
   //cout << "term vel is "<<(1.0/6.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,1)/(VISCOSITY*DENS)<<endl;
 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<endl;

   Nmisc::boundaryCylinder(ps,CYLINDER_ORIGIN,CYLINDER_RADIUS,CYLINDER_HEIGHT);
   vect origin = CYLINDER_ORIGIN;
   origin[2] += 1.5*PSEP;
   Nmisc::sphCylinder(ps,origin);

   p.tag = ps.size()+1;
   p.r = CYLINDER_ORIGIN;
   p.r[2] = CYLINDER_HEIGHT;
   p.dens = DEM_DENS;
   p.h = DEM_RADIUS/2;
   p.mass = (4.0/3.0)*PI*pow(DEM_RADIUS,3)*DEM_DENS;
   p.v = 0.0;
   p.vhat = p.v;
   p.iam = dem;
   ps.push_back(p);

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
      ioFile.setFilename(filename.c_str(),&globals);
      ioFile.writeGlobals(0,&globals);
      ioFile.writeRestart(0,vps[i],&globals);
      globals.mpiRank = 0;
   }
   cout << "Writing Global data to file..."<<endl;
   //ioFile.writeGlobals(0,&globals);
   ioFile.writeDomain(0,vprocDomain,vprocNeighbrs);

}

