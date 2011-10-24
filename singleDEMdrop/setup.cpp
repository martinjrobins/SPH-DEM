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

   //do boundaries here

   vect points[4];
   points[0] = BMIN[0],BMAX[1];
   points[1] = BMIN[0],BMIN[1];
   points[2] = BMAX[0],BMIN[1];
   points[3] = BMAX[0],BMAX[1];

   vect normals[4];
   normals[0] = 1,0;
   normals[1] = 0,1;
   normals[2] = -1,0;
   normals[3] = -1,0;

   bool concave[4];
   concave[0] = true;
   concave[1] = true;
   concave[2] = true;
   concave[3] = true;
   
   //Nmisc::boundaryLine(ps,points,normals,concave,4,boundary);
   int numBound1,numBound2;
   if (PERIODIC[0]) {
      numBound1 = 1;
      numBound2 = 0;
   } else {
      numBound1 = 3;
      numBound2 = 3;
   }

   for (int i=1-numBound1;i<=NX+numBound2-1;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=1-numBound1;j<=NY+numBound2-1;j++) {
         if ((j<=0)&&(i<=0)&&(!PERIODIC[0])&&(!PERIODIC[1])) {
            p.iam = sphBoundary;
            p.norm1[0] = -i*PSEP;
            p.norm1[1] = -j*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[0] = 1.0/sqrt(2);
               p.norm1[1] = 1.0/sqrt(2);
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((j<=0)&&(i>=NX)&&(!PERIODIC[0])&&(!PERIODIC[1])) {
            p.iam = sphBoundary;
            p.norm1[0] = (NX-i)*PSEP;
            p.norm1[1] = -j*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[0] = -1.0/sqrt(2);
               p.norm1[1] = 1.0/sqrt(2);
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((j>=NY)&&(i<=0)&&(!PERIODIC[0])&&(!PERIODIC[1])) {
            p.iam = sphBoundary;
            p.norm1[0] = -i*PSEP;
            p.norm1[1] = (NY-j)*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[0] = 1.0/sqrt(2);
               p.norm1[1] = -1.0/sqrt(2);
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((j>=NY)&&(i>=NX)&&(!PERIODIC[0])&&(!PERIODIC[1])) {
            p.iam = sphBoundary;
            p.norm1[0] = (NX-i)*PSEP;
            p.norm1[1] = (NY-j)*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[0] = -1.0/sqrt(2);
               p.norm1[1] = -1.0/sqrt(2);
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((j<=0)&&(!PERIODIC[0])) {
            p.iam = sphBoundary;
            p.norm1[0] = 0;
            p.norm1[1] = -j*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[1] = 1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((i<=0)&&(!PERIODIC[1])) {
            p.iam = sphBoundary;
            p.norm1[0] = -i*PSEP;
            p.norm1[1] = 0;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[0] = 1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((i>=NX)&&(!PERIODIC[1])) {
            p.iam = sphBoundary;
            p.norm1[0] = (NX-i)*PSEP;
            p.norm1[1] = 0;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[0] = -1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else {
            p.iam = sph;
         }
         p.tag = ps.size()+1;
         p.r = (i)*PSEP+BMIN[0],(j)*PSEP+BMIN[1];
         p.dens = DENS;
         p.mass = PSEP*PSEP*DENS;
         //cout <<p.mass<<endl;
         p.h = H;
         p.v = 0.0;
         p.vhat = p.v;
         p.alpha = ALPHA;
         if ((p.iam==sph)&&(p.r[1]>BMIN[1]+(BMAX[1]-BMIN[1])*0.97)) continue;
         ps.push_back(p);
      }
   }
   p.tag = ps.size()+1;
   p.r[0] = (BMAX[0]-BMIN[0])/2;
   p.r[1] = 1.00*(BMAX[1]-BMIN[1]);
   p.dens = DEM_DENS;
   p.h = DEM_RADIUS/2;
//#ifdef _2D_
   //p.mass = PI*pow(DEM_RADIUS,2)*DEM_DENS;
//#else
   p.mass = (4.0/3.0)*PI*pow(DEM_RADIUS,3)*DEM_DENS;
   //cout << p.mass<<endl;
//#endif
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

