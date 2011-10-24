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

void setNormals(vect &norm1,vect &norm3,vect newNorm1,vect newNorm2,int &count) {
   if (count>0) {
      vect normtmp = (norm1*count + newNorm1)/(count+1);
      norm3 = cross(norm1,newNorm1);
      norm1 = normtmp;
   } else {
      norm1 = newNorm1;
      norm3 = newNorm2;
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

   for (int i=0;i<NX;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=0;j<=NY;j++) {
         for (int k=0;k<=NZ;k++) {
            p.tag = ps.size()+1;
            p.r = (i)*PSEP+RMIN[0],(j)*PSEP+RMIN[1],k*PSEP+RMIN[2];
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
            if (k<=0) {
               p.iam = boundary;
               newNorm1 = 0.0,0.0,1.0;
               newNorm2 = 0.0,1.0,0.0;
               setNormals(p.norm1,p.norm3,newNorm1,newNorm2,count);
            } 
            if (i<=0) {
               p.iam = boundary;
               newNorm1 = 1,0,0;
               newNorm2 = 0,0,1;
               setNormals(p.norm1,p.norm3,newNorm1,newNorm2,count);
            } 
            if (j <= 0) {
               p.iam = boundary;
               newNorm1 = 0,1,0;
               newNorm2 = 1,0,0;
               setNormals(p.norm1,p.norm3,newNorm1,newNorm2,count);
            } 
            if (j >= NY) {
               p.iam = boundary;
               newNorm1 = 0,-1,0;
               newNorm2 = 1,0,0;
               setNormals(p.norm1,p.norm3,newNorm1,newNorm2,count);
            }
            if (p.iam==boundary) {
               p.norm1 = p.norm1/len(p.norm1);
               p.norm3 = p.norm3/len(p.norm3);
            }
            if ((p.iam==sph)&&(p.r[2]>H1+H2)&&(p.r[0]>L1)) continue;
            if ((p.iam==sph)&&(p.r[2]>H1)&&(p.r[0]>L1+L2)) continue;
            ps.push_back(p);
         }
      }
   }

   for (int i=NX;i<=NX+ceil(L4/PSEP);i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=0;j<=NY;j++) {
         for (int k=0;k<=NZ;k++) {
            p.tag = ps.size()+1;
            p.r = (i)*PSEP+RMIN[0],(j)*PSEP+RMIN[1],k*PSEP+RMIN[2];
            p.dens = DENS;
            p.mass = pow(PSEP,NDIM)*DENS;
            p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
            p.v = 0.0;
            p.vhat = p.v;
            p.iam = sph;
            p.norm1 = 0;
            p.norm3 = 0;
            if (k==0) {
               p.iam = boundary;
               p.norm1 = 0,0,1;
               p.norm3 = 0,1,0;
               p.norm2 = 0;
            } 
            if ((p.iam==sph)&&((p.r[1] < PSEP + 0.5*(p.r[0]-NX*PSEP)/L4)||(p.r[1] > 1-PSEP-0.5*(p.r[0]-NX*PSEP)/L4))) continue;
            if ((p.iam==sph)&&(p.r[2] > H1)) continue;
            ps.push_back(p);
         }
      }
   }

   vect norm1;
   norm1 = L4,0.5,0;
   vect norm2;
   norm2 = 0,0,H1+H2+H3+PSEP;
   vect initPos;
   initPos = PSEP*NX,0,0;
   vect neighbrNorm1,neighbrNorm2,neighbrNorm3,neighbrNorm4;
   neighbrNorm1 = 0,1,0;
   neighbrNorm2 = 0;
   vect tmp;
   tmp = L4,-0.5,0;
   vect tmpNorm = tmp/len(tmp);
   vect norm2Norm = norm2/len(norm2);
   neighbrNorm3 = cross(tmpNorm,norm2Norm);
   neighbrNorm4 = 0.0;

   Nmisc::boundaryPlane(ps,norm1,norm2,neighbrNorm1,neighbrNorm2,neighbrNorm3,neighbrNorm4,initPos,PSEP,true);

   norm1 = tmp;
   initPos = PSEP*NX,RMAX[1],0;
   neighbrNorm1 = 0,-1,0;
   neighbrNorm2 = 0,0,0;
   neighbrNorm3 = 0.0;
   neighbrNorm4 = 0.0;
   Nmisc::boundaryPlane(ps,norm1,norm2,neighbrNorm1,neighbrNorm2,neighbrNorm3,neighbrNorm4,initPos,PSEP,false);
             
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
   //ioFile.writeGlobals(0,&globals);
   //ioFile.writeDomain(0,vprocDomain,vprocNeighbrs);

}

