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
#include <gsl/gsl_chebyshev.h>
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

   int numBound1,numBound2;
   if (PERIODIC[0]||GHOST[0]) {
      numBound1 = 1;
      numBound2 = 0;
   } else {
      numBound1 = 4;
      numBound2 = 4;
   }

      
   for (int i=1-numBound1;i<=NX+numBound2-1;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=1-numBound1;j<=NY+numBound2-1;j++) {
         p.tag = ps.size()+1;
         if (PERIODIC[0]||GHOST[0]) {
            p.r = (i+0.5)*PSEP+RMIN[0],(j+0.5)*PSEP+RMIN[1];
         } else {
            p.r = i*PSEP+RMIN[0],j*PSEP+RMIN[1];
         }
         p.dens = DENS;
         p.mass = PSEP*PSEP*DENS;
         p.h = H;
         //p.v = -p.r[1]*BOX_GLOB_ANGVEL,p.r[0]*BOX_GLOB_ANGVEL;
         const double r1 = sqrt(pow(p.r[0]-CENTRE1[0],2)+pow(p.r[1]-CENTRE1[1],2));
         const double r2 = sqrt(pow(p.r[0]-CENTRE2[0],2)+pow(p.r[1]-CENTRE2[1],2));
         p.v[0] = -0.5*abs(OMEGA_E1)*(p.r[1]-CENTRE1[1])*exp(-pow(r1/R0,2))+0.5*abs(OMEGA_E1)*(p.r[1]-CENTRE2[1])*exp(-pow(r2/R0,2));
         p.v[1] = -0.5*abs(OMEGA_E1)*(p.r[0]-CENTRE1[0])*exp(-pow(r1/R0,2))+0.5*abs(OMEGA_E1)*(p.r[0]-CENTRE2[0])*exp(-pow(r2/R0,2));
         p.vhat = p.v;

         if ((j<=0)&&(i<=0)&&(!PERIODIC[0])&&(!PERIODIC[1])&&(!GHOST[0])) {
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
         } else if ((j<=0)&&(i>=NX)&&(!PERIODIC[0])&&(!PERIODIC[1]&&(!GHOST[0]))) {
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
         } else if ((j>=NY)&&(i<=0)&&(!PERIODIC[0])&&(!PERIODIC[1]&&(!GHOST[0]))) {
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
         } else if ((j>=NY)&&(i>=NX)&&(!PERIODIC[0])&&(!PERIODIC[1])&&(!GHOST[0])) {
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
         } else if ((j<=0)&&(!PERIODIC[0])&&(!GHOST[0])) {
            p.iam = sphBoundary;
            p.norm1[0] = 0;
            p.norm1[1] = -j*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[1] = 1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((j>=NY)&&(!PERIODIC[0])&&(!GHOST[0])) {
            p.iam = sphBoundary;
            p.norm1[0] = 0;
            p.norm1[1] = (NY-j)*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[1] = -1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((i<=0)&&(!PERIODIC[1])&&(!GHOST[0])) {
            p.iam = sphBoundary;
            p.norm1[0] = -i*PSEP;
            p.norm1[1] = 0;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[0] = 1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((i>=NX)&&(!PERIODIC[1])&&(!GHOST[0])) {
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
         ps.push_back(p);
      }
   }


   double totKEfinal=0;
   double rmsVfinal=0;
   double maxVLenfinal = 0;
#ifdef FORCING
   maxVLen = 1;
   totKE = 1;
#endif
   for (int i=0;i<ps.size();i++) {
      totKEfinal += 0.5*ps[i].mass*len2(ps[i].v);
      rmsVfinal += len2(ps[i].v);
      if (len(ps[i].v)>maxVLenfinal) maxVLenfinal = len(ps[i].v);
   }
   rmsVfinal /= ps.size();
   cout <<" Final total KE = "<<totKEfinal<<endl;
   cout <<" Final max velocity = "<<maxVLenfinal<<endl;
   cout <<" Final rms velocity = "<<rmsVfinal<<endl;
   //maxVLen *= sqrt(4.0/sumMV2);
   //cout << "Final max V len = "<<maxVLen<<endl;

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

