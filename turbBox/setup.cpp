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
   cout <<"damp time = "<<BOX_GLOB_TIME<<endl;
 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<endl;

   const double beta = 100;

   Array<double,NDIM> XChebCoeff(CHEB_N+1,CHEB_N+1);
   Array<double,NDIM> YChebCoeff(CHEB_N+1,CHEB_N+1);

   gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd1);
   
   //gsl_rng_set(rng,1);
   gsl_rng_set(rng,time(NULL));

   for (int i=0;i<=CHEB_N;i++) {
      for (int j=0;j<=CHEB_N;j++) {
         double sigma;
         if ((i==CHEB_N)||(j==CHEB_N)) {
            sigma = 0.0;
         } else {
            sigma = (i/double(1+pow(i/8.0,4)))*(j/double(1+pow(j/8.0,4)));
         }
         XChebCoeff(i,j) = gsl_ran_gaussian(rng,sigma);
         YChebCoeff(i,j) = gsl_ran_gaussian(rng,sigma);
      }
   }
   gsl_rng_free(rng);

   double maxVLen = 0.0;
   double sumMV2 = 0.0;

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
      for (int j=1-numBound1;j<=NY+numBound2-1;j++) {
         p.tag = ps.size()+1;
         if (PERIODIC[0]) {
            p.r = (i+0.5)*PSEP+RMIN[0],(j+0.5)*PSEP+RMIN[1];
         } else {
            p.r = i*PSEP+RMIN[0],j*PSEP+RMIN[1];
         }
         p.dens = DENS;
         p.mass = PSEP*PSEP*DENS;
         p.h = H;
         //p.v = -p.r[1]*BOX_GLOB_ANGVEL,p.r[0]*BOX_GLOB_ANGVEL;
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
         } else if ((j>=NY)&&(!PERIODIC[0])) {
            p.iam = sphBoundary;
            p.norm1[0] = 0;
            p.norm1[1] = (NY-j)*PSEP;
            p.dist = len(p.norm1);
            if (p.dist==0) {
               p.norm1[1] = -1.0;
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
            //p.v[0] += gsl_cheb_eval(csX,p.r[0]);
            //p.v[1] += gsl_cheb_eval(csY,p.r[1]);
#ifndef FORCING
            vect vRand = 0.0;
            for (int a=0;a<=CHEB_N;a++) {
               for (int b=0;b<=CHEB_N;b++) {
                  vRand[0] += XChebCoeff(a,b)*cos(a*acos(p.r[0]))*cos(b*acos(p.r[1]));
                  vRand[1] += YChebCoeff(a,b)*cos(a*acos(p.r[0]))*cos(b*acos(p.r[1]));
               }
            }
            vRand *= (1-exp(-beta*pow(1-pow(p.r[0],2),2)));
            vRand *= (1-exp(-beta*pow(1-pow(p.r[1],2),2)));
            if (len(vRand)>maxVLen) maxVLen = len(vRand);
            sumMV2 += p.mass*len2(vRand);
            p.v += vRand;
#endif
         }
         p.alpha = ALPHA;
         ps.push_back(p);
      }
   }

   //setup arrays
   const double DX = (RMAX[0]-RMIN[0])/NX;
   const double DY = (RMAX[1]-RMIN[1])/NY;
   //const double DX = PSEP;
   //const double DY = PSEP;

   Array<double,NDIM> s(NX+1,NY+1);
   Array<double,NDIM> phi(NX+1,NY+1);
   Array<Cparticle *,NDIM> pp(NX+1,NY+1);

   int pn = 0;
   for (int i=-3;i<=NX+3;i++) {
      for (int j=-3;j<=NY+3;j++) {
         if ((j>=0)&&(j<=NY)&&(i>=0)&&(i<=NX)) {
               pp(i,j) = &(ps[pn]);
         }
         pn++;
      }
   }

   for (int i=0;i<=NX;i++) {
      for (int j=0;j<=NY;j++) {
         phi(i,j) = 0.0;
         if ((j==0)||(j==NY)||(i==0)||(i==NX)) {
            s(i,j) = 0;
         } else {
            s(i,j) = (pp(i-1,j)->v[0]-pp(i+1,j)->v[0])/(2*DX) + (pp(i,j-1)->v[1]-pp(i,j+1)->v[1])/(2*DY);
         }
      }
   }

   cout << "maximum -div v = "<<max(abs(s))<<"maximum velocity magnitude = "<<maxVLen<<endl;

   cout << "solving elliptical equation..."<<endl;
#ifdef FORCING
   for (int iteration=0;iteration<0;iteration++) {
#else
   for (int iteration=0;iteration<60000;iteration++) {
#endif
      double maxDelPhi = 0;
      cout << "iteration "<<iteration<<endl;
      for (int i=1;i<NX;i++) {
         for (int j=1;j<NY;j++) {
            double oldPhi = phi(i,j);
            if ((i==1)&&(j==1)) {
               //phi(i,j) = (phi(i+1,j)+phi(i,j+1))/2;
               //phi(i,j) = 0;
               //phi(i,j) = phi(i+1,j);
               phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j+1);
            } else if ((i==1)&&(j==NY-1)) {
               //phi(i,j) = (phi(i+1,j)+phi(i,j-1))/2;
               //phi(i,j) = 0;
               //phi(i,j) = phi(i+1,j);
               phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j-1);
            } else if ((i==NX-1)&&(j==1)) {
               //phi(i,j) = (phi(i-1,j)+phi(i,j+1))/2;
               //phi(i,j) = 0;
               //phi(i,j) = phi(i,j+1);
               phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j+1);
            } else if ((i==NX-1)&&(j==NY-1)) {
               //phi(i,j) = (phi(i-1,j)+phi(i,j-1))/2;
               //phi(i,j) = 0;
               //phi(i,j) = phi(i,j-1);
               phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j-1);
            } else if (i==1) {
               //phi(i,j) = (phi(i+1,j)+phi(i,j-1)+phi(i,j+1))/3;
               //phi(i,j) = 0;
               //phi(i,j) = -pow(DX,2)*s(i,j) + phi(i+1,j);
               phi(i,j) = -(pow(DY,2)*pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(2*pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*(phi(i,j+1)+phi(i,j-1));
            } else if (i==NX-1) {
               //phi(i,j) = (phi(i-1,j)+phi(i,j-1)+phi(i,j+1))/3;
               //phi(i,j) = 0;
               //phi(i,j) = -pow(DX,2)*s(i,j) + phi(i-1,j);
               phi(i,j) = -(pow(DY,2)*pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(2*pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*(phi(i,j+1)+phi(i,j-1));
            } else if (j==1) {
               //phi(i,j) = (phi(i-1,j)+phi(i+1,j)+phi(i,j+1))/3;
               //phi(i,j) = 0;
               //phi(i,j) = -pow(DY,2)*s(i,j) + phi(i,j+1);
               phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+2*pow(DY,2)))*(phi(i+1,j)+phi(i-1,j))+(pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*phi(i,j+1);
            } else if (j==NY-1) {
               //phi(i,j) = (phi(i-1,j)+phi(i+1,j)+phi(i,j-1))/3;
               //phi(i,j) = 0;
               //phi(i,j) = -pow(DY,2)*s(i,j) + phi(i,j-1);
               phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+2*pow(DY,2)))*(phi(i+1,j)+phi(i-1,j))+(pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*phi(i,j-1);
            } else {
               phi(i,j) = -0.5*(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+0.5*(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*(phi(i+1,j)+phi(i-1,j))+0.5*(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*(phi(i,j+1)+phi(i,j-1));
            }
            double delPhi = phi(i,j)-oldPhi;
            if (delPhi > maxDelPhi) maxDelPhi = delPhi;
            pp(i,j)->vort = phi(i,j);
         }
      }
      cout << "\tmaxDelPhi = "<<maxDelPhi<<endl;
      if (maxDelPhi==0) break;
   }

   for (int i=1;i<NX;i++) {
      phi(i,0) = phi(i,1);
      phi(i,NY) = phi(i,NY-1);
   }
   for (int j=1;j<NY;j++) {
      phi(0,j) = phi(1,j);
      phi(NX,j) = phi(NX-1,j);
   }
   phi(0,0) = phi(1,1);
   phi(0,NY) = phi(1,NY-1);
   phi(NX,0) = phi(NX-1,1);
   phi(NX,NY) = phi(NX-1,NY-1);

   for (int i=1;i<NX;i++) {
      for (int j=1;j<NY;j++) {
         vect u;
         u[0] = (phi(i+1,j)-phi(i-1,j))/(2*DX);
         u[1] = (phi(i,j+1)-phi(i,j-1))/(2*DY);
         pp(i,j)->alpha = s(i,j)-((phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/pow(DX,2) + (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/pow(DY,2));
         pp(i,j)->v += u;
         pp(i,j)->v *= (1-exp(-beta*pow(1-pow(pp(i,j)->r[0],2),2)));
         pp(i,j)->v *= (1-exp(-beta*pow(1-pow(pp(i,j)->r[1],2),2)));
      }
   }
   maxVLen = 0;
   double totKE = 0;
   for (int i=0;i<=NX;i++) {
      for (int j=0;j<=NY;j++) {
         if ((j==0)||(j==NY)||(i==0)||(i==NX)) {
            s(i,j) = 0;
         } else {
            s(i,j) = (pp(i-1,j)->v[0]-pp(i+1,j)->v[0])/(2*DX) + (pp(i,j-1)->v[1]-pp(i,j+1)->v[1])/(2*DY);
         }
         pp(i,j)->alpha = s(i,j);
         if (len(pp(i,j)->v)>maxVLen) maxVLen = len(pp(i,j)->v);
	 totKE += 0.5*pp(i,j)->mass*len2(pp(i,j)->v);
      }
   }

   cout << "maximum -div v = "<<max(abs(s))<<"maximum velocity magnitude = "<<maxVLen<<endl;

   double totKEfinal=0;
   double rmsVfinal=0;
   double maxVLenfinal = 0;
#ifdef FORCING
   maxVLen = 1;
   totKE = 1;
#endif
   for (int i=0;i<ps.size();i++) {
#ifdef SPIN
      ps[i].v *= MAX_VRAND/maxVLen;
      if (ps[i].iam==sphBoundary) {
         vect br = ps[i].r+ps[i].norm1*ps[i].dist;
         ps[i].v[0] += -br[1]*BOX_GLOB_ANGVEL;
         ps[i].v[1] += br[0]*BOX_GLOB_ANGVEL;
         ps[i].vhat[0] += -ps[i].r[1]*BOX_GLOB_ANGVEL;
         ps[i].vhat[1] += ps[i].r[0]*BOX_GLOB_ANGVEL;
      } else {
         ps[i].v[0] += -ps[i].r[1]*BOX_GLOB_ANGVEL;
         ps[i].v[1] += ps[i].r[0]*BOX_GLOB_ANGVEL;
         ps[i].vhat = ps[i].v;
      }
#else
      ps[i].v *= sqrt(KE_REF/totKE);
      ps[i].vhat = ps[i].v;
#endif
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

