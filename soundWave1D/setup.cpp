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
   cout << "Creating length with sides: rmax = ["<<RMAX[0]<<"] rmin = ["<<RMIN[0]<<"]"<<endl;
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

   int numBound1,numBound2;
   if (PERIODIC[0]) {
      numBound1 = 1;
      numBound2 = 0;
   } else {
      numBound1 = 4;
      numBound2 = 4;
   }

   for (int i=1-numBound1;i<=NX+numBound2-1;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
         p.tag = ps.size()+1;
         if (PERIODIC[0]) {
            p.r[0] = (i+0.5)*PSEP+RMIN[0];
         } else {
            p.r[0] = i*PSEP+RMIN[0];
         }
         p.dens = DENS;
         p.mass = PSEP*DENS;
         p.mass += gsl_ran_gaussian(rng,0.2*p.mass);
         p.h = H;
         p.v = 0,0;
         //p.v = -p.r[1]*BOX_GLOB_ANGVEL,p.r[0]*BOX_GLOB_ANGVEL;
         //const vect middle = (RMAX[0]-RMIN[0])/2.0 + RMIN[0];
         //const vect size = (RMAX[0]-RMIN[0])/7.0;
         //const double sigma = (RMAX[0]-RMIN[0])/14.0;
         //if (all(p.r>=middle)&&all(p.r<=middle+size)) {
         //   p.v = sin(2.0*PI*(1.0/size)*(p.r-middle));
         //} else {
         //   p.v = 0;
         //}
         //if (all(p.r>RMIN[0])&&all(p.r<RMAX[0])) {
         //   p.v = 0.1*exp(-len2(p.r-middle)/(2.0*pow(sigma,2)));
            //p.v = 0.1*sin(2.0*PI*3.0*p.r);
         //} else {
         //   p.v = 0.0;
         //}

         if ((i<=0)&&(!PERIODIC[0])) {
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
            p.v[0] = 0.1*sin(PI*p.r[0]/PSEP);
            p.v[0] += 0.1*sin(10.0*PI*p.r[0]/(RMAX[0]-RMIN[0]));
            cout << "p.v  = "<<p.v<<endl;
         }
         p.alpha = ALPHA;
         ps.push_back(p);
   }


   cout << "Total number of particles = " << ps.size() << endl;

   CglobalVars globals;


   vector<vector<double> > vprocDomain(globals.mpiSize);
   vector<Array<int,NDIM> > vprocNeighbrs(globals.mpiSize);
   vector<particleContainer > vps;
   vectInt split;
   split = 1;
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


   //Write restart file for John
   /*ofstream fo("restartJohn.dat");
   fo <<"h = "<<H<<" mass = "<<ps[0].mass<<" dens = "<<ps[0].dens<<" psep = "<<PSEP<<" Re = "<<REYNOLDS_NUMBER<<" kinematic viscosity = "<<VISCOSITY<<" n = "<<ps.size()<<endl;
   fo <<"r_x r_y v_x v_y flag(0=fluid,1=boundary)"<<endl;
   for (int i=0;i<ps.size();i++) {
      fo << ps[i].r[0]<<' '<<ps[i].r[1]<<' '<<ps[i].v[0]<<' '<<ps[i].v[1]<<' '<<ps[i].iam<<endl;
   }*/

}

