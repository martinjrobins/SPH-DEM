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

   gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd1);
   gsl_rng_set(rng,1);

   int numBound1,numBound2,numBound3,numBound4;
   if (PERIODIC[0]) {
      numBound1 = 1;
      numBound2 = 0;
   } else {
      numBound1 = 4;
      numBound2 = 4;
   }
   if (PERIODIC[1]) {
      numBound3 = 1;
      numBound4 = 0;
   } else {
      numBound3 = 4;
      numBound4 = 4;
   }

      
   vect middle,size;
   for (int i=1;i<NDIM;i++) {
      middle[i] = (RMAX[i]-RMIN[i])/2.0 + RMIN[i];
      size[i] = (RMAX[i]-RMIN[i])/7.0;
   }
   //const double sigma = (RMAX[0]-RMIN[0])/14.0;
   const double sigma = 0.5*PSEP;
   for (int i=1-numBound1;i<=NX+numBound2-1;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
         if (i==20) continue;
      for (int j=1-numBound3;j<=NY+numBound4-1;j++) {
         p.tag = ps.size()+1;
         
         if (i==21) {
            if (PERIODIC[0]) {
               p.r[0] = i*PSEP+RMIN[0];
            } else {
               p.r[0] = (i-0.5)*PSEP+RMIN[0];
            }
            p.mass = 2*PSEP*PSEP*DENS;
         } else {
            if (PERIODIC[0]) {
               p.r[0] = (i+0.5)*PSEP+RMIN[0];
            } else {
               p.r[0] = i*PSEP+RMIN[0];
            }
            p.mass = PSEP*PSEP*DENS;
         }
         if (PERIODIC[1]) {
            p.r[1] = (j+0.5)*PSEP+RMIN[1];
         } else {
            p.r[1] = j*PSEP+RMIN[1];
         }
         p.dens = DENS;
         //p.mass += gsl_ran_gaussian(rng,0.2*p.mass);

         p.h = H;
         
         p.v = 0,0;

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
         } else if ((j<=0)&&(!PERIODIC[1])) {
            p.iam = sphBoundary;
            p.norm1[0] = 0;
            p.norm1[1] = -j*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[1] = 1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((j>=NY)&&(!PERIODIC[1])) {
            p.iam = sphBoundary;
            p.norm1[0] = 0;
            p.norm1[1] = (NY-j)*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[1] = -1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((i<=0)&&(!PERIODIC[0])) {
            p.iam = sphBoundary;
            p.norm1[0] = -i*PSEP;
            p.norm1[1] = 0;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[0] = 1.0;
            } else {
               p.norm1 = p.norm1/p.dist;
            }
         } else if ((i>=NX)&&(!PERIODIC[0])) {
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

            //p.v = 0.1*exp(-len2(p.r-middle)/(2.0*pow(sigma,2)));
            //p.v[0] = 0.1*exp(-pow(p.r[0]-middle[0],2)/(2.0*pow(sigma,2)));
            p.v[0] = 0.1*cos(PI*i);
            cout <<"p.v = "<<p.v<<" i = "<<i<<" j = "<<j<<" sigma = "<<sigma<<" p.r = "<<p.r<<" middle = "<<middle<<" cos(PI*i) = "<<cos(PI*i)<<endl;
            //p.v = 0.1*sin(2.0*PI*3.0*p.r);
            //p.v[0] = 0.1*sin(2.0*PI*3.0*p.r[0]);

            //vect oneOne;
            //oneOne[0] = 1.0/sqrt(2.0);
            //oneOne[1] = 1.0/sqrt(2.0);
            //vect oneNegOne;
            //oneNegOne[0] = -1.0/sqrt(2.0);
            //oneNegOne[1] = 1.0/sqrt(2.0);
            //vect rmov = p.r;
            //rmov[0] -= RMIN[0];
            //rmov[1] -= RMIN[1];
            //const double xhat = len(rmov-dot(rmov,oneOne)*oneOne);
            //p.v = 0.1*(1.0/sqrt(2.0))*oneNegOne*(exp(-pow(xhat,2)/(2.0*pow(sigma,2))) + exp(-pow(xhat-sqrt(2),2)/(2.0*pow(sigma,2))));
            
            
            //p.v[0] += gsl_cheb_eval(csX,p.r[0]);
            //p.v[1] += gsl_cheb_eval(csY,p.r[1]);
         }
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

