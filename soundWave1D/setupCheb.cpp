#include "customConstants.h"
#include "particle.h"
#include "io_data_vtk.h"
#include "misc.h"
#include "sph.h"
#include "sphIncompress.h"
#include "customOutput.h"
#include "dataLL.h"
#include <vector>
#include <iostream>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define CHEB_N 65
#define CHEB_SEED 0

int main() {
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

   const double beta = 100;

   Array<double,NDIM> XChebCoeff(CHEB_N+1,CHEB_N+1);
   Array<double,NDIM> YChebCoeff(CHEB_N+1,CHEB_N+1);

   gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd1);
   gsl_rng_set(rng,CHEB_SEED);

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

   for (int i=-2;i<=NX+2;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=-2;j<=NY+2;j++) {
         p.tag = ps.size()+1;
         p.r = i*PSEP+RMIN[0],j*PSEP+RMIN[1];
         p.dens = DENS;
         p.mass = PSEP*PSEP*DENS;
         p.h = H;
         p.v = -p.r[1]*BOX_GLOB_ANGVEL,p.r[0]*BOX_GLOB_ANGVEL;

         /*if ((j==0)&&(i==0)) {
            p.iam = boundary;
            p.norm = 1.0,1.0;
            p.norm = p.norm/len(p.norm);
         } else if ((j==0)&&(i==NX)) {
            p.iam = boundary;
            p.norm = -1.0,1.0;
            p.norm = p.norm/len(p.norm);
         } else if ((j==NY)&&(i==0)) {
            p.iam = boundary;
            p.norm = 1.0,-1.0;
            p.norm = p.norm/len(p.norm);
         } else if ((j==NY)&&(i==NX)) {
            p.iam = boundary;
            p.norm = -1.0,-1.0;
            p.norm = p.norm/len(p.norm);
         } else if (j==0) {
            p.iam = boundary;
            p.norm = 0.0,1.0;
         } else if (j==NY) {
            p.iam = boundary;
            p.norm = 0.0,-1.0;
         } else if (i==0) {
            p.iam = boundary;
            p.norm = 1.0,0.0;
         } else if (i==NX) {
            p.iam = boundary;
            p.norm = -1.0,0.0;
         } else {
            p.norm = 0,0;
            p.iam = sph;
         }
         */
         if ((j<=0)||(j>=NY)||(i<=0)||(i>=NX)) {
            p.iam = sphBoundary;
         } else {
            p.iam = sph;
            //p.v[0] += gsl_cheb_eval(csX,p.r[0]);
            //p.v[1] += gsl_cheb_eval(csY,p.r[1]);
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
   for (int i=-2;i<=NX+2;i++) {
      for (int j=-2;j<=NY+2;j++) {
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
   for (int iteration=0;iteration<40000;iteration++){
      double maxDelPhi = 0;
      cout << "iteration "<<iteration<<endl;
      for (int i=1;i<NX;i++) {
         for (int j=1;j<NY;j++) {
            double oldPhi = phi(i,j);
            if ((i==1)&&(j==1)) {
               //phi(i,j) = (phi(i+1,j)+phi(i,j+1))/2;
               //phi(i,j) = 0;
               phi(i,j) = phi(i+1,j);
               //phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j+1);
            } else if ((i==1)&&(j==NY-1)) {
               //phi(i,j) = (phi(i+1,j)+phi(i,j-1))/2;
               //phi(i,j) = 0;
               phi(i,j) = phi(i+1,j);
               //phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j-1);
            } else if ((i==NX-1)&&(j==1)) {
               //phi(i,j) = (phi(i-1,j)+phi(i,j+1))/2;
               //phi(i,j) = 0;
               phi(i,j) = phi(i,j+1);
               //phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j+1);
            } else if ((i==NX-1)&&(j==NY-1)) {
               //phi(i,j) = (phi(i-1,j)+phi(i,j-1))/2;
               //phi(i,j) = 0;
               phi(i,j) = phi(i,j-1);
               //phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(pow(DX,2)+pow(DY,2)))*phi(i,j-1);
            } else if (i==1) {
               //phi(i,j) = (phi(i+1,j)+phi(i,j-1)+phi(i,j+1))/3;
               //phi(i,j) = 0;
               phi(i,j) = -pow(DX,2)*s(i,j) + phi(i+1,j);
               //phi(i,j) = -(pow(DY,2)*pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(2*pow(DX,2)+pow(DY,2)))*phi(i+1,j)+(pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*(phi(i,j+1)+phi(i,j-1));
            } else if (i==NX-1) {
               //phi(i,j) = (phi(i-1,j)+phi(i,j-1)+phi(i,j+1))/3;
               //phi(i,j) = 0;
               phi(i,j) = -pow(DX,2)*s(i,j) + phi(i-1,j);
               //phi(i,j) = -(pow(DY,2)*pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*s(i,j)+(pow(DY,2)/(2*pow(DX,2)+pow(DY,2)))*phi(i-1,j)+(pow(DX,2)/(2*pow(DX,2)+pow(DY,2)))*(phi(i,j+1)+phi(i,j-1));
            } else if (j==1) {
               //phi(i,j) = (phi(i-1,j)+phi(i+1,j)+phi(i,j+1))/3;
               //phi(i,j) = 0;
               phi(i,j) = -pow(DY,2)*s(i,j) + phi(i,j+1);
               //phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+2*pow(DY,2)))*(phi(i+1,j)+phi(i-1,j))+(pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*phi(i,j+1);
            } else if (j==NY-1) {
               //phi(i,j) = (phi(i-1,j)+phi(i+1,j)+phi(i,j-1))/3;
               //phi(i,j) = 0;
               phi(i,j) = -pow(DY,2)*s(i,j) + phi(i,j-1);
               //phi(i,j) = -(pow(DY,2)*pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*s(i,j)+(pow(DY,2)/(pow(DX,2)+2*pow(DY,2)))*(phi(i+1,j)+phi(i-1,j))+(pow(DX,2)/(pow(DX,2)+2*pow(DY,2)))*phi(i,j-1);
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
   for (int i=0;i<=NX;i++) {
      for (int j=0;j<=NY;j++) {
         if ((j==0)||(j==NY)||(i==0)||(i==NX)) {
            s(i,j) = 0;
         } else {
            s(i,j) = (pp(i-1,j)->v[0]-pp(i+1,j)->v[0])/(2*DX) + (pp(i,j-1)->v[1]-pp(i,j+1)->v[1])/(2*DY);
         }
         if (len(pp(i,j)->v)>maxVLen) maxVLen = len(pp(i,j)->v);
      }
   }

   cout << "maximum -div v = "<<max(abs(s))<<"maximum velocity magnitude = "<<maxVLen<<endl;

   for (int i=0;i<ps.size();i++) {
      //ps[i].v *= sqrt(4.0/sumMV2);
      ps[i].v *= MAX_VRAND/maxVLen;
   }
   //maxVLen *= sqrt(4.0/sumMV2);
   //cout << "Final max V len = "<<maxVLen<<endl;

   cout << "Total number of particles = " << ps.size() << endl;

   CglobalVars globals;

   CdataLL *data = new CdataLL(ps,globals);
   CsphIncompress sph(data);
   CcustomOutput customOutput(data,globals.time);

   cout << "Calculating Output stuff.."<<endl;
   sph.calcOutputVars();
   customOutput.calcOutput(globals.time);

   cout << "Opening files for writing..."<<endl;
   Cio_data_vtk ioFile("initData");
   cout << "Writing Restart data to file..."<<endl;
	ioFile.writeRestart(0,ps);
   cout << "Writing Global data to file..."<<endl;
   ioFile.writeGlobals(0,&globals);


   //Write restart file for John
   ofstream fo("restartJohn.dat");
   fo <<"h = "<<H<<" mass = "<<ps[0].mass<<" dens = "<<ps[0].dens<<" psep = "<<PSEP<<" Re = "<<REYNOLDS_NUMBER<<" kinematic viscosity = "<<VISCOSITY<<endl;
   fo <<"r_x r_y v_x v_y flag(0=fluid,1=boundary)"<<endl;
   for (int i=0;i<ps.size();i++) {
      fo << ps[i].r[0]<<' '<<ps[i].r[1]<<' '<<ps[i].v[0]<<' '<<ps[i].v[1]<<' '<<ps[i].iam<<endl;
   }



}

