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


void demPorusCircle(vector<Cparticle> &ps,const vect origin,const double radius,const double psep) {
   Cparticle p;
   const double NR = ceil(radius/psep);
   const double RPSEP = radius/NR;
   const vect n(1,0,0);
   for (int i=0;i<=NR;i++) {
      const double r = i*RPSEP;
      double nc = ceil(2*PI*r/psep);
      if (nc == 0) nc++;
      const double TPSEP = 2.0*PI/nc;
      for (int j=0;j<nc;j++) {
         const double theta = j*TPSEP;
         vect newN;
         newN[0] = n[0]*cos(theta)-n[1]*sin(theta);
         newN[1] = n[0]*sin(theta)+n[1]*cos(theta);
         newN[2] = 0;
         p.r = r*newN + origin;
         p.tag = ps.size()+1;
         p.dens = DEM_DENS;
         p.mass = (4.0/3.0)*PI*pow(DEM_RADIUS,3)*DEM_DENS;
         p.h = H;
         p.v = 0,0,0;
         p.iam = demBoundary;
         ps.push_back(p);
      }
   }
}


void demPorousCylinder(vector<Cparticle> &ps,const vect origin,const double radius, const double height, const double discPorosity) {
   const double demVol = (4.0/3.0)*PI*pow(DEM_RADIUS,3);
   const double diskVol = PI*pow(radius,2)*height;
   const double solidVol = diskVol*(1-discPorosity);
   const double numDem = solidVol/demVol;
   cout <<" number of dem particles needed = "<<numDem<<endl;
   const double numDemPerUnitVol = numDem/diskVol;
   cout <<" number of dem particle per unit vol = "<<numDemPerUnitVol<<endl;
   const double avePsep = pow(numDemPerUnitVol,-1.0/NDIM);
   cout <<" calculated avePsep for porus disks is "<<avePsep<<". This is "<<avePsep/(DEM_RADIUS*2.0)<<" times the dem diameter"<<endl;

   const int psSize = ps.size();
   const int NH = static_cast<int>(height/avePsep);
   const double hpsep = height/NH;
   vect newOrigin = origin;
   for (int i=0;i<NH;i++) {
      newOrigin[2] = origin[2] + i*hpsep;
      demPorusCircle(ps,newOrigin,radius,avePsep);
   }
   const int n = ps.size()-psSize;
   const double finalPorosity = (diskVol-n*demVol)/diskVol;
   cout <<"After putting in this porus disk, the final porosity is "<<finalPorosity<<". The goal porosity was "<<discPorosity<<endl;
}

void myBoundaryCircle(vector<Cparticle> &ps,const vect &origin,const double radiusMin,const double radiusMax,enum iamTypes iam) {
   Cparticle p;
   const double dRadius = radiusMax-radiusMin;
   const double NR = ceil(dRadius/(BFAC*PSEP));
   double RPSEP = dRadius/NR;
   if (dRadius==0.0) RPSEP = 0;
   const vect n(1,0,0);
   for (int i=0;i<=NR;i++) {
      const double r = i*RPSEP+radiusMin;
      double nc = int(ceil(2.0*PI*r/(BFAC*PSEP)));
      if (nc == 0) nc++;
      const double TPSEP = 2.0*PI/nc;
      for (int j=0;j<nc;j++) {
         const double theta = j*TPSEP;
         vect newN;
         newN[0] = n[0]*cos(theta)-n[1]*sin(theta);
         newN[1] = n[0]*sin(theta)+n[1]*cos(theta);
         newN[2] = 0;
         p.r = r*newN + origin;
         p.tag = ps.size()+1;
         p.dens = DENS;
         p.mass = pow(BFAC*PSEP,NDIM)*DENS;
         p.h = H;
         p.v = 0,0,0;
         p.vhat = p.v;
         p.iam = iam;
         p.norm3 = newN;
         if (i==NR) {
            p.norm1 = -cos(theta),-sin(theta),1;
            p.norm1 = p.norm1/len(p.norm1);
            //p.norm2[0] = -cos(theta);
            //p.norm2[1] = -sin(theta);
            //p.norm2[2] = 0;
            p.norm2 = 0,0,0;
            p.norm3[2] += 1;
            p.norm3 = p.norm3/len(p.norm3);
         } else {
            p.norm1 = 0,0,1;
            p.norm2 = 0,0,0;
         }
         ps.push_back(p);
      }
   }
}



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
   cout << "this run will need "<<MAX_NUM_PARTICLES_PER_CPU*sizeof(Cparticle)<<" bytes"<<endl;
   cout <<"inlet vel is "<<INFLOW_VEL<<endl;

   cout << "PSEP = "<<PSEP<<endl;
   cout << "alpha = "<<ALPHA<<endl;
   cout << "viscosity = "<<VISCOSITY<<endl;
   cout << "maxtime = "<<MAXTIME<<endl;
   cout << "term vel (3d) is "<<(2.0/9.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,2)/(VISCOSITY*DENS)<<endl;
   cout << "term vel (2d) is "<<(1.0/3.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,2)/(VISCOSITY*DENS)<<endl;
   //cout << "term vel is "<<(1.0/6.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,1)/(VISCOSITY*DENS)<<endl;

   cout << "gamma must be less than" <<2.0*sqrt(DEM_K*DEM_MIN_REDUCED_MASS)<<endl;
#ifdef LUBRICATION
   cout << "k must be greater than" <<pow(0.5*LUB_GAMMA,2)/DEM_MIN_REDUCED_MASS<<endl;
#else
   cout << "k must be greater than" <<pow(0.5*DEM_GAMMA,2)/DEM_MIN_REDUCED_MASS<<endl;
#endif
   cout << "Konmr ="<<Konmr<<endl;
   cout << "gammaOnmr = "<<gammaOnmr<<endl;
   cout << "K = "<<DEM_K<<endl;
#ifdef LUBRICATION
   cout << "GAMMA  = "<<LUB_GAMMA<<endl;
#else
   cout << "GAMMA  = "<<DEM_GAMMA<<endl;
#endif
   

   cout <<"THETS ="<<THETS<<endl;
   cout <<"DEM_TIMESTEP="<<DEM_TIMESTEP<<endl;
   cout <<"CR = "<<CR<<endl;

 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   double tdem = Nsph::demCondition();
#ifdef LUBRICATION
   cout <<"coeff of restitution = "<<exp(-tdem*50*LUB_GAMMA/(2.0*DEM_MIN_REDUCED_MASS))<<endl;
#else
   cout <<"coeff of restitution = "<<exp(-tdem*50*DEM_GAMMA/(2.0*DEM_MIN_REDUCED_MASS))<<endl;
#endif
   double liqtdem = Nsph::liqDemCondition();
   cout <<"dem ts = "<<tdem<<endl;
   cout <<"liq dem ts = "<<liqtdem<<endl;
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<int((MAXTIME/tdem)+1)<<" steps according to the DEM condition"<<int((MAXTIME/liqtdem)+1)<<" steps according to the LIQ DEM condition"<<endl;

   cout <<"reduced mass = "<<DEM_MIN_REDUCED_MASS<<endl;
cout <<" tdem = "<< (1.0/50.0)*PI*sqrt(DEM_MIN_REDUCED_MASS)/sqrt(DEM_K-pow(0.5*DEM_GAMMA,2)/DEM_MIN_REDUCED_MASS) << endl;
   cout <<"particle reynolds number = "<<2.0*DEM_RADIUS*INFLOW_VEL/VISCOSITY<<endl;
   vect normal = 0.0;
   normal[2] = 1.0;
   vect newOrigin = CYLINDER_ORIGIN;
#ifdef SPH_INLET
   Nmisc::boundaryCircle(ps,CYLINDER_ORIGIN,INLET_RADIUS+3*PSEP,CYLINDER_RADIUS,normal,0,1);
#else
   Nmisc::boundaryCircle(ps,CYLINDER_ORIGIN,INLET_RADIUS,CYLINDER_RADIUS,normal);
#endif
   
   CglobalVars globals;
   globals.procNeighbrs = -1;
   for (int j=0;j<NDIM;j++) {
      globals.procDomain[2*j] = RMIN[j]-0.25*(RMAX[j]-RMIN[j]);
      globals.procDomain[2*j+1] = RMAX[j]+0.25*(RMAX[j]-RMIN[j]);
   }
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &(globals.mpiSize));
   MPI_Comm_rank(MPI_COMM_WORLD, &(globals.mpiRank));
   ps.reserve(MAX_NUM_PARTICLES_PER_CPU);
   cout <<"creating new data structure"<<endl;
   CdataLL *data = new CdataLL(ps,globals,true);
   
#ifdef LIQ_DEM
   cout <<"adding dem particles"<<endl;
   
   cout <<"going to add "<<NDEM<<" dem particles at random locations...."<<endl;
   newOrigin[2] = CYLINDER_ORIGIN[2]+1.0*DEM_RADIUS;
   cout <<"adding cylinder at origin"<<newOrigin<<endl;
   Nmisc::addRandomDEMparticlesCyl(data,newOrigin,CYLINDER_RADIUS-1.0*DEM_RADIUS,CYLINDER_HEIGHT-2.0*DEM_RADIUS,NDEM);
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
      if (p->iam == dem) {
         p->mass = INIT_DEM_MASS;
         p->dens = INIT_DEM_DENS;
      }
   }
#endif
   //newOrigin = CYLINDER_ORIGIN;
   //origin[2] += 1.6*PSEP;
   //Nmisc::sphCylinder(ps,origin);
   
   
   cout <<"adding cylinder at origin"<<newOrigin<<endl;
   newOrigin = CYLINDER_ORIGIN;
   newOrigin[2] += BFAC*PSEP;
   const double newHeight = CYLINDER_HEIGHT*(1+TOP_BUFFER) - BFAC*PSEP;
   Nmisc::boundaryCylinderNoTopBottom(ps,newOrigin,CYLINDER_RADIUS,newHeight);
   newOrigin[2] = CYLINDER_ORIGIN[2]-INLET_HEIGHT;
   cout <<"adding cylinder at origin"<<newOrigin<<endl;
#ifdef SPH_INLET
   Nmisc::sphBoundaryCylinderNoTopBottom(ps,newOrigin,INLET_RADIUS,INLET_HEIGHT+PSEP);
#else
   Nmisc::boundaryCylinderNoTopBottom(ps,newOrigin,INLET_RADIUS,INLET_HEIGHT);
#endif
   newOrigin[2] = CYLINDER_ORIGIN[2]+CYLINDER_HEIGHT*(1+TOP_BUFFER);
   normal[2] = -1.0;
   cout <<"adding cylinder at origin"<<newOrigin<<endl;
#ifdef SPH_INLET
   Nmisc::boundaryCircle(ps,newOrigin,INLET_RADIUS+3*PSEP,CYLINDER_RADIUS,normal,0,1);
#else
   Nmisc::boundaryCircle(ps,newOrigin,INLET_RADIUS,CYLINDER_RADIUS,normal);
#endif
   cout <<"adding cylinder at origin"<<newOrigin<<endl;
#ifdef SPH_INLET
   Nmisc::sphBoundaryCylinderNoTopBottom(ps,newOrigin,INLET_RADIUS,INLET_HEIGHT);
#else
   newOrigin[2] += BFAC*PSEP;
   Nmisc::boundaryCylinderNoTopBottom(ps,newOrigin,INLET_RADIUS,INLET_HEIGHT);
#endif
   //demPorousCylinder(ps,CYLINDER_ORIGIN,CYLINDER_RADIUS,DEM_DISK_HEIGHT,DEM_DISK_POROSITY);
   //newOrigin[2] = CYLINDER_ORIGIN[2]+CYLINDER_HEIGHT-DEM_DISK_HEIGHT;
   //demPorousCylinder(ps,newOrigin,CYLINDER_RADIUS,DEM_DISK_HEIGHT,0.5);


#ifdef WET_START
   const double totalDEMVol = NDEM*(4.0/3.0)*PI*pow(DEM_RADIUS,3);
   //const double newLiqDens = DENS;
   const double liqVol = PI*pow(CYLINDER_RADIUS,2)*CYLINDER_HEIGHT*0.90;
   const double newLiqDens = DENS*(liqVol-totalDEMVol)/(liqVol);
   const double pmass = pow(PSEP,NDIM)*DENS;
   const double newPSEP = pow(pmass/newLiqDens,1.0/NDIM);
   const double newNX = (BMAX[0]-BMIN[0])/newPSEP;
   const double newNY = newNX;
   const double newNZ = CYLINDER_HEIGHT*(1+TOP_BUFFER)/newPSEP;
   //const double newLiqDens = DENS;
   cout << "porosity = "<<1.0-totalDEMVol/(liqVol)<<endl;
   cout << "after adding dem particles, new liquid density is "<<newLiqDens<<endl; 

   cout <<"adding liquid particles"<<endl;
   for (int i=0;i<=newNX;i++) {
         cout << "\rParticle ("<<i<<","<<"0"<<"). Generation "<<((i+2)*(NY+4))/double((NX+4)*(NY+4))*100<<"\% complete"<<flush;
      for (int j=0;j<=newNY;j++) {
         for (int k=0;k<=newNZ;k++) {
            p.tag = ps.size()+1;
            p.r = (i)*newPSEP+BMIN[0],j*newPSEP+BMIN[1],(k+1)*newPSEP+CYLINDER_ORIGIN[2];
            if ((p.r[2]<CYLINDER_ORIGIN[2]+0.5*PSEP)||(p.r[2]>CYLINDER_ORIGIN[2]+CYLINDER_HEIGHT*(1+TOP_BUFFER)-0.5*PSEP)||(p.r[0]*p.r[0]+p.r[1]*p.r[1]>pow(CYLINDER_RADIUS-0.5*PSEP,2))) continue;
            p.dens = newLiqDens;
            p.mass = pmass;
            p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
            p.v = 0.0;
            p.vhat = p.v;
            p.iam = sph;
            ps.push_back(p);
         }
      }
   }

   for (int k=0;k<int(INLET_HEIGHT/PSEP)+1;k++) {
      vect tmp = CYLINDER_ORIGIN;
      tmp[2] -= INLET_HEIGHT-k*PSEP;
      if (k<4) {
         myBoundaryCircle(ps,tmp,0.5*(NIX%2)*PSEP,NIX*PSEP/2.0,sphBoundary);
      } else {
         myBoundaryCircle(ps,tmp,0.5*(NIX%2)*PSEP,NIX*PSEP/2.0,sph);
      }

      tmp[2] = CYLINDER_ORIGIN[2]+CYLINDER_HEIGHT*(1+TOP_BUFFER)+k*PSEP;
      myBoundaryCircle(ps,tmp,0.5*(NIX%2)*PSEP,NIX*PSEP/2.0,sph);
   }

#endif
   /*
   p.tag = ps.size()+1;
   p.r = CYLINDER_ORIGIN;
   p.r[2] = CYLINDER_HEIGHT;
   p.dens = DEM_DENS;
   p.h = DEM_RADIUS/2;
   p.mass = (4.0/3.0)*PI*pow(DEM_RADIUS,3)*DEM_DENS;
   p.v = 0.0;
   p.vhat = p.v;
   p.iam = dem;
   //ps.push_back(p);
   */

   cout << "Total number of particles = " << ps.size() << endl;

   //CglobalVars globals;

   vector<vector<double> > vprocDomain(globals.mpiSize);
   vector<Array<int,NDIM> > vprocNeighbrs(globals.mpiSize);
   vector<particleContainer > vps;
   vectInt split;
   split = NCPU_X,NCPU_Y,NCPU_Z;
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

