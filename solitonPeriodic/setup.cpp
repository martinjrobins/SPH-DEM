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

void setNormals(Cparticle &p,vect newNorm1,vect newNorm2,int &count) {
   if (count==0) {
      p.norm1 = newNorm1;
      p.norm3 = newNorm2;
   //} else if (count==1) {
   //   p.norm2 = newNorm1;
   //   p.norm3 = cross(p.norm1,p.norm2);
   //}
   } else {
      p.norm3 = cross(p.norm1,newNorm1);
      p.norm1 = (p.norm1*count + newNorm1)/(count+1);
   }
   count++;
} 


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
   cout << "speed of sound = "<<SPSOUND<<endl;
   cout << "prb = "<<PRB<<endl;
   cout << "number of particles on side = "<<NX<<endl;
   cout << "max num partilces = "<<MAX_NUM_PARTICLES_PER_CPU<<endl;

   cout << "PSEP = "<<PSEP<<endl;
   cout << "alpha = "<<ALPHA<<endl;
   cout << "viscosity = "<<VISCOSITY<<endl;
   cout << "maxtime = "<<MAXTIME<<endl;
 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<endl;

   int NXwall = int((L1+L2)/PSEP);
   cout <<"ENDPERIODICCPU = "<<ENDPERIODICCPU<<endl;

   for (int i=0;i<NXP;i++) {
      cout <<"i = "<<i<<"p.r[0] = "<<(i+0.5)*PSEP+RMIN[0]<<endl;
      for (int j=(NY/2.0)-PY;j<(NY/2.0)+PY;j++) {
         for (int k=0;k<NZ;k++) {
            p.tag = STARTTAG;
            p.r = (i+0.5)*PSEP+RMIN[0],(j+0.5)*PSEP+RMIN[1],(k+0.5)*PSEP+RMIN[2];
            p.dens = DENS;
            p.mass = pow(PSEP,NDIM)*DENS;
            p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
            p.v = 0.0;
            p.vhat = p.v;
            p.norm2 = 0;
            p.norm1 = 0;
            p.norm3 = 0;
#ifdef FLUID_GATE
            if ((i <= NXwall)&&(i>=NXwall-gateWidth)) {
               p.iam = sphBoundary;
#else
            if ((i <= NXwall)&&(i>=NXwall-gateWidth)) {
            //if ((i == NXwall)) {
               p.iam = boundary;
               if (i==NXwall) {
                  p.norm1[0] = 1;
               } else if (i==NXwall-gateWidth) {
                  p.norm1[0] = -1;
               } else if (k>0) {
                  continue;
               }
               p.norm3[1] = 1;
               if (k==0) {
                  if ((i==NXwall)||(i==NXwall-gateWidth)) {
                     p.norm2[2] = -1;
                  } else {
                     p.norm1[0] = 0;
                     p.norm1[2] = -1;
                  }
               }
#endif
            } else {
               p.iam = sph;
            }
            if ((p.iam==sph)&&(p.r[2]>H1+H2)&&(p.r[0]>L1)) continue;
            if ((p.iam==sph)&&(p.r[2]>H1)&&(p.r[0]>L1+L2)) continue;
            ps.push_back(p);
         }
      }
   }

   cout << "Finished first section" << endl;
   cout << "Total number of particles = " << ps.size() << endl;

   for (int i=NXP;i<NX;i++) {
      cout <<"i = "<<i<<"p.r[0] = "<<(i+0.5)*PSEP+RMIN[0]<<endl;
      for (int j=0;j<NY;j++) {
         for (int k=0;k<NZ;k++) {
            p.tag = ENDTAG;
            p.r = (i+0.5)*PSEP+RMIN[0],(j+0.5)*PSEP+RMIN[1],(k+0.5)*PSEP+RMIN[2];
            p.dens = DENS;
            p.mass = pow(PSEP,NDIM)*DENS;
            p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
            p.v = 0.0;
            p.vhat = p.v;
            int count = 0;
            p.norm2 = 0;
            vect newNorm1,newNorm2;
            p.iam = sph;
            p.norm1 = 0;
            p.norm3 = 0;
            if ((p.iam==sph)&&(p.r[2]>H1+H2)&&(p.r[0]>L1)) continue;
            if ((p.iam==sph)&&(p.r[2]>H1)&&(p.r[0]>L1+L2)) continue;
            ps.push_back(p);
         }
      }
   }

   cout << "Total number of particles = " << ps.size() << endl;
   cout <<"creating tip"<<endl;
   for (int i=NX;i<=NX+ceil(L5/PSEP);i++) {
      cout <<"i = "<<i<<"p.r[0] = "<<(i+0.5)*PSEP+RMIN[0]<<endl;
      for (int j=0;j<=NY;j++) {
         for (int k=0;k<=NZ;k++) {
            p.tag = ENDTAG;
            p.r = (i+0.5)*PSEP+RMIN[0],(j+0.5)*PSEP+RMIN[1],(k+0.5)*PSEP+RMIN[2];
            p.dens = DENS;
            p.mass = pow(PSEP,NDIM)*DENS;
            p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
            p.v = 0.0;
            p.vhat = p.v;
            p.iam = sph;
            p.norm1 = 0;
            p.norm3 = 0;
            if ((p.iam==sph)&&((p.r[1] < PSEP + 0.5*WIDTH*(p.r[0]-NX*PSEP)/L5)||(p.r[1] > WIDTH-PSEP-0.5*WIDTH*(p.r[0]-NX*PSEP)/L5))) continue;
            if ((p.iam==sph)&&(p.r[2] > H1)) continue;
            ps.push_back(p);
         }
      }
   }

   cout << "Total number of particles = " << ps.size() << endl;
   vect norm1;
   norm1 = L5,WIDTH/2.0,0;
   vect norm1Norm = norm1/len(norm1);
   //norm1 = norm1 - PSEP*norm1/len(norm1);
   vect norm2;
   norm2 = 0,0,WALLH+PSEP;
   vect initPos;
   initPos = PSEP*NX,0,-0.5*PSEP;
   //initPos = initPos + norm1*PSEP/len(norm1);
   vect neighbrNorm1,neighbrNorm2,neighbrNorm3,neighbrNorm4;
   neighbrNorm1 = 0,0,0;
   neighbrNorm2 = 0;
   vect tmp;
   tmp = L5,-WIDTH/2.0,0;
   vect tmpNorm = tmp/len(tmp);
   vect norm2Norm = norm2/len(norm2);
   neighbrNorm3 = cross(tmpNorm,norm2Norm);
   neighbrNorm4 = 0.0;

   Nmisc::boundaryPlane(ps,norm1,norm2,neighbrNorm1,false,neighbrNorm2,false,neighbrNorm3,true,neighbrNorm4,false,initPos,PSEP,true);

   norm1 = tmp;
   initPos = PSEP*NX,RMAX[1]-0.0001*PSEP,-0.5*PSEP;
   neighbrNorm1 = cross(norm2Norm,norm1Norm);
   neighbrNorm2 = 0,0,0;
   neighbrNorm3 = 0.0;
   neighbrNorm4 = 0.0;
   Nmisc::boundaryPlane(ps,norm1,norm2,neighbrNorm1,false,neighbrNorm2,false,neighbrNorm3,false,neighbrNorm4,false,initPos,PSEP,false);

             
   cout << "Total number of particles = " << ps.size() << endl;


   CglobalVars globals;
   vector<vector<double> > vprocDomain(globals.mpiSize);
   vector<Array<int,NDIM> > vprocNeighbrs(globals.mpiSize);
   vector<particleContainer > vps;
   vectInt split;
   split = NCPU,1,1;
   particleContainer pps;
   for (int i=0;i<ps.size();i++) {
      pps.push_back(ps[i]);
   }
   Nmisc::splitDomain(pps,split,vps,vprocDomain,vprocNeighbrs);

   const double width1 = (NXP)*PSEP/(ENDPERIODICCPU+1);
   const double width2 = (L4+L5)/(NCPU-ENDPERIODICCPU-1);

   for (int i=0;i<=ENDPERIODICCPU;i++) {
      vprocDomain[i][0] = RMIN[0]+i*width1;
      vprocDomain[i][1] = RMIN[0]+(i+1)*width1;
      cout <<"cpu = "<<i<<" max boundary = "<<vprocDomain[i][1]<<endl;
   }
   for (int i=ENDPERIODICCPU+1;i<NCPU;i++) {
      vprocDomain[i][0] = RMIN[0]+(ENDPERIODICCPU+1)*width1 + (i-ENDPERIODICCPU-1)*width2;
      vprocDomain[i][1] = RMIN[0]+(ENDPERIODICCPU+1)*width1 + (i-ENDPERIODICCPU)*width2;
      cout <<"cpu = "<<i<<" max boundary = "<<vprocDomain[i][1]<<endl;
   }
   vprocDomain[NCPU-1][1] = RMAX[0]+0.25*(RMAX[0]-RMIN[0]);

   for (int i=0;i<NCPU;i++) {
      vps[i].clear();
   }
   cout <<"dividing particles..."<<endl;
   Nmisc::divideParticlesToCPUs(pps,vps,vprocDomain);
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
   //ioFile.writeGlobals(0,&globals);
   //ioFile.writeDomain(0,vprocDomain,vprocNeighbrs);

}

