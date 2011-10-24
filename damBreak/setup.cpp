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

   vect points[5];
   points[0] = BMIN[0],BMAX[1];
   points[1] = BMIN[0],BMIN[1];
   points[2] = BMAX[0],BMIN[1];
   points[3] = BMAX[0],BMAX[1];
   points[4] = BMIN[0],BMAX[1];

   vect normals[5];
   normals[0] = 1,0;
   normals[1] = 0,1;
   normals[2] = -1,0;
   normals[3] = 0,-1;
   normals[4] = 1,0;

   bool concave[5];
   concave[0] = true;
   concave[1] = true;
   concave[2] = true;
   concave[3] = true;
   concave[4] = true;
   
   Nmisc::boundaryLine(ps,points,normals,concave,5,boundary);


   cout << "Total number of boundary particles = " << ps.size() << endl;

   //do fluid here
   for (int i=2;i<=NX;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=2;j<=NY;j++) {

         p.r = i*PSEP+BMIN[0],j*PSEP+BMIN[1];

         p.tag = ps.size()+1;
         p.dens = DENS;
         p.mass = PSEP*PSEP*DENS;
         p.h = H;
         p.v = 0,0;
         p.iam = sph;
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

