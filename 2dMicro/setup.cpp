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


int main(int argc, char *argv[]) {
   if (argc != 2) {
      cout << "Usage: setup outfilename" << endl;
      return(-1);
   }
   string filename = argv[1];

   vector<Cparticle> ps;
   Cparticle p;
   cout << "Creating box with sides: rmax = ["<<BMAX[0]<<" "<<BMAX[1]<<"] rmin = ["<<BMIN[0]<<" "<<BMIN[1]<<"]"<<endl;
   cout << "Particle Separation = "<<PSEP<<endl;
   cout << "Reynolds Number = "<<REYNOLDS_NUMBER<<endl;
   cout << "Density = "<<DENS<<endl;
   cout << "number of particles on side = "<<NX<<endl;
   cout << "alpha = "<<ALPHA<<endl;
   cout << "viscosity = "<<VISCOSITY<<endl;
   cout << "maxtime = "<<MAXTIME<<endl;
 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<endl;

   //do boundaries here

   vect points2[4];
   double cx = (BMAX[0]-BMIN[0])/2.0;
   points2[0] = tan(THETA)*PSEP+cx-WW1/2.0,PSEP;
   points2[1] = tan(THETA)*H1+cx-WW1/2.0,H1;
   points2[2] = cx+WW1/2.0-H1*tan(THETA),H1;
   points2[3] = cx+WW1/2.0,0;

   vect normals2[4];
   normals2[0] = -cos(THETA),sin(THETA);
   normals2[1] = 0,1;
   normals2[2] = cos(THETA),sin(THETA);
   normals2[3] = -cos(THETA),sin(THETA);
   
   bool concave2[4];
   concave2[0] = false;
   concave2[1] = false;
   concave2[2] = false;
   concave2[3] = false;
   
   Nmisc::boundaryLine(ps,points2,normals2,concave2,4,boundaryobj1);



   vect points[9];
   points[0] = H1*tan(THETA),0;
   points[1] = BMAX[0]-H1*tan(THETA),0;
   points[2] = BMAX[0],H1;
   points[3] = BMAX[0],BMAX[1]-H2;
   points[4] = BMAX[0]-H2*tan(THETA),BMAX[1];
   points[5] = H2*tan(THETA),BMAX[1];
   points[6] = 0,BMAX[1]-H2;
   points[7] = 0,H1;
   points[8] = H1*tan(THETA),0;

   vect normals[9];
   normals[0] = 0,1;
   normals[1] = -cos(THETA),sin(THETA);
   normals[2] = -1,0;
   normals[3] = -cos(THETA),-sin(THETA);
   normals[4] = 0,-1;
   normals[5] = cos(THETA),-sin(THETA);
   normals[6] = 1,0;
   normals[7] = cos(THETA),sin(THETA);
   normals[8] = 0,1;

   bool concave[9];
   concave[0] = true;
   concave[1] = true;
   concave[2] = true;
   concave[3] = true;
   concave[4] = true;
   concave[5] = true;
   concave[6] = true;
   concave[7] = true;
   concave[8] = true;

   Nmisc::boundaryLine(ps,points,normals,concave,9,boundary);

   #ifndef NO_TOP
   vect points3[4];
   points3[0] = PSEP*tan(THETA)+cx-WW2/2.0,BMAX[1]-PSEP;
   points3[1] = H2*tan(THETA)+cx-WW2/2.0,BMAX[1]-H2;
   points3[2] = cx+WW2/2.0-H2*tan(THETA),BMAX[1]-H2;
   points3[3] = cx+WW2/2.0,BMAX[1];

   vect normals3[4];
   normals3[0] = -cos(THETA),-sin(THETA);
   normals3[1] = 0,-1;
   normals3[2] = cos(THETA),-sin(THETA);
   normals3[3] = -cos(THETA),-sin(THETA);
   
   bool concave3[4];
   concave3[0] = false;
   concave3[1] = false;
   concave3[2] = false;
   concave3[3] = false;
   
   Nmisc::boundaryLine(ps,points3,normals3,concave3,4,boundaryobj2);
#endif


   cout << "Total number of boundary particles = " << ps.size() << endl;

   //do fluid here
   for (int i=1;i<NX;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=1;j<NY;j++) {

         p.r = i*PSEP+BMIN[0],j*PSEP+BMIN[1];

         bool outside = false;
         double rtanT = 1.0/tan(THETA);
         const double buffer = 0.9;
         outside |= p.r[1]<=rtanT*p.r[0] + H1-BMAX[0]*rtanT + buffer*PSEP;
         outside |= p.r[1]>=-rtanT*p.r[0] + BMAX[1]-H2+BMAX[0]*rtanT - buffer*PSEP;
         outside |= p.r[1]>=rtanT*p.r[0] + BMAX[1]-H2 - buffer*PSEP;
         outside |= p.r[1]<=-rtanT*p.r[0] + H1 + buffer*PSEP;
         outside |= (p.r[1]<=H1+buffer*PSEP)&&(p.r[0]>=cx-WW1/2.0-buffer*PSEP)&&(p.r[0]<=cx+WW1/2.0+buffer*PSEP)&&(p.r[1]<=rtanT*(p.r[0]-cx+WW1/2.0)+buffer*PSEP)&&(p.r[1]<=-rtanT*(p.r[0]-cx-WW1/2.0)+buffer*PSEP);
#ifndef NO_TOP
         outside |= (p.r[1]>=BMAX[1]-H2-buffer*PSEP)&&(p.r[0]>=cx-WW2/2.0-buffer*PSEP)&&(p.r[0]<=cx+WW2/2.0+buffer*PSEP)&&(p.r[1]>=-rtanT*(p.r[0]-cx+WW2/2.0)+BMAX[1]-buffer*PSEP)&&(p.r[1]>=rtanT*(p.r[0]-cx-WW2/2.0)+BMAX[1]-buffer*PSEP);
#endif

         if (!outside) {
            p.tag = ps.size()+1;
            p.dens = DENS;
            p.mass = PSEP*PSEP*DENS;
            p.h = H;
            p.v = 0,0;
            p.iam = sph;
            p.alpha = ALPHA;
            ps.push_back(p);
         }
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

