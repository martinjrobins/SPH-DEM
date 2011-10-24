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
   cout << "Reynolds Number = "<<REYNOLDS_NUMBER<<endl;
   cout << "Density = "<<DENS<<endl;
   cout << "speed of sound = "<<SPSOUND<<endl;
   cout << "prb = "<<PRB<<endl;
   cout << "number of particles on side = "<<NX<<endl;
   cout << "max num partilces = "<<MAX_NUM_PARTICLES_PER_CPU<<endl;
   cout << "sph boundary density is "<<SPHBOUNDARYDENS<<endl;
   cout << "VMAX is "<<VMAX<<endl;

   cout << "PSEP = "<<PSEP<<endl;
   cout << "alpha = "<<ALPHA<<endl;
   cout << "viscosity = "<<VISCOSITY<<endl;
   cout << "maxtime = "<<MAXTIME<<endl;
   cout << "term vel (3d) is "<<VREF<<endl;
   cout << "term vel (stokes) is "<<(2.0/9.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,2)/(VISCOSITY*DENS)<<endl;
   cout << "DEM_RADIUS = "<<DEM_RADIUS<<endl;
   cout << "Other dens = "<<D2<<endl;
   cout << "dens drop= "<<DENS_DROP[0]<<endl;
 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   double tdem = Nsph::demCondition();
   cout <<"dem ts = "<<tdem<<endl;
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<int((MAXTIME/tdem)+1)<<" steps according to the DEM condition"<<endl;

   for (int i=0;i<NX;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=0;j<NY;j++) {
         for (int k=0;k<NZ;k++) {
            p.tag = ps.size()+1;
            p.r = (i+0.5)*PSEP+RMIN[0],(j+0.5)*PSEP+RMIN[1],(k+0.5)*PSEP+RMIN[2];
            if (p.tag==11) cout <<p.r<<endl;
            p.dens = DENS;
            p.mass = pow(PSEP,NDIM)*DENS;
            p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
            p.v = 0.0;
            p.vhat = p.v;
            ps.push_back(p);
         }
      }
   }
   int numSPH = ps.size();

   CglobalVars globals;
   globals.procNeighbrs = 0;
   globals.procNeighbrs(1,1,1) = -1;
   for (int j=0;j<NDIM;j++) {
      if (!PERIODIC[j]) {
         globals.procDomain[2*j] = RMIN[j]-0.25*(RMAX[j]-RMIN[j]);
      } else {
         globals.procDomain[2*j] = RMIN[j];
      }
      if (!PERIODIC[j]) {
         globals.procDomain[2*j+1] = RMAX[j]+0.25*(RMAX[j]-RMIN[j]);
      } else {
         globals.procDomain[2*j+1] = RMAX[j];
      }
   }
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &(globals.mpiSize));
   MPI_Comm_rank(MPI_COMM_WORLD, &(globals.mpiRank));
   ps.reserve(MAX_NUM_PARTICLES_PER_CPU);
   cout <<"creating new data structure"<<endl;
   CdataLL *data = new CdataLL(ps,globals,true);
   cout <<"adding dem particles"<<endl;
#ifdef MANY_PARTICLES
   double min[3];
   double max[3];
   for (int i=0;i<3;i++) {
      if (i==2) {
         min[i] = RMIN[i];
         max[i] = RMAX[i];
      } else {
         min[i] = RMIN[i];
         max[i] = RMAX[i];
      }
   }
#ifdef GRID_OF_DEM
   Nmisc::addGridDEMparticles(ps,min,max,POROSITY);
#else
   Nmisc::addRandomDEMparticles(data,min,max,POROSITY);
#endif
#else
   p.tag = ps.size()+1;
   p.r[0] = (RMAX[0]-RMIN[0])/2;
   p.r[1] = (RMAX[1]-RMIN[1])/2;
   p.r[2] = RMAX[2]-4.0*PSEP;
   p.dens = DEM_DENS;
   p.h = DEM_RADIUS;
   p.mass = DEM_VOL*DEM_DENS;
   p.v = 0.0;
   p.vhat = p.v;
   p.iam = dem;
   ps.push_back(p);
#endif



   const double totalDEMVol = (ps.size()-numSPH)*(4.0/3.0)*PI*pow(DEM_RADIUS,3);
   //const double newLiqDens = DENS;
   const double newLiqDens = DENS*((RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0])*(RMAX[2]-RMIN[2])-totalDEMVol)/((RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0])*(RMAX[2]-RMIN[2]));
   cout << "porosity = "<<1.0-totalDEMVol/((RMAX[0]-RMIN[0])*(RMAX[1]-RMIN[1])*(RMAX[2]-RMIN[2]))<<endl;
   cout << "after adding dem particles, new liquid density is "<<newLiqDens<<endl; 

   for (vector<Cparticle>::iterator i=ps.begin();i!=ps.end();i++) {
      if (i->iam==sph) {
         i->dens = newLiqDens;
         i->mass = pow(PSEP,NDIM)*newLiqDens;
      }
   }


   cout << "Total number of particles = " << ps.size() << endl;


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

