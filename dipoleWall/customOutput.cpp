#include "customOutput.h"
#include "sph.h"
#include "misc.impl.h"
#include "dataLL.impl.h"

using namespace Nsph;

inline void calcEnstrophyAndAveDens(Cparticle &p,CglobalVars &g) {
   g.custom[0] += 0.5*p.mass*len2(p.vort)/p.dens;
   g.custom[1] += p.mass*dot(p.vhat,p.fv)/p.dens;
   g.custom[2] += p.mass*dot(p.vhat,p.ff)/p.dens;
   g.custom[3] += p.mass*dot(p.vhat,p.f)/p.dens;
   g.custom[4] += p.mass*p.pdr2*p.dddt/p.dens;
   g.custom[5] += p.mass*dot(p.vhat,p.fp)/p.dens;
   g.custom[6] += 0.5*g.dt*p.mass*dot(p.ff,p.ff)/p.dens;
   g.custom[7] += p.mass*dot(p.fv,p.fv)/p.dens;
   g.custom[8] += p.mass*p.press/p.dens;
}


void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {

#ifdef VORT_LEASTSQUARES
   data->neighboursGroup<calcVortLeastSquares,ifSph>();
#else
   data->neighboursGroup<calcVortSPHSum,ifSph>();
#endif
   data->globals.custom[0] = 0;
   data->globals.custom[1] = 0;
   data->globals.custom[2] = 0;
   data->globals.custom[3] = 0;
   data->globals.custom[4] = 0;
   data->globals.custom[5] = 0;
   data->globals.custom[6] = 0;
   data->globals.custom[7] = 0;
   data->globals.custom[8] = 0;
   data->traverse<calcEnstrophyAndAveDens,ifSphOrSphBoundary>();
   //data->globals.custom[1] /= data->numParticles();
}

inline bool ifInteriorSph(Cparticle &p) {
   return (p.iam == sph)&&(abs(p.r[0])<RMAX[0]*0.75)&&(abs(p.r[1])<RMAX[1]*0.75);
}

inline void calcMaxVandDens(Cparticle &p,CglobalVars &g) {
   double v = len(p.v);
   if (v > g.maxV) g.maxV = v;
   g.aveDens += p.dens;
}

inline void calcTheVarDens(Cparticle &p,CglobalVars &g) {
   g.varDens += pow(p.dens-g.aveDens,2);
}

inline double velSpec(Cparticle &p,CglobalVars &g) {
   if (PERIODIC[0]) {
      return (0.5*p.mass*len2(p.v));
   } else {
      return (0.5*p.mass*len2(p.v)-(g.eKE/g.nSph));
   }
}

#ifdef AVE_VELOCITY
inline double aveVSpec(Cparticle &p,CglobalVars &g) {
   if (PERIODIC[0]) {
      return (0.5*p.mass*len2(p.aveV));
   } else {
      return (0.5*p.mass*len2(p.aveV)-(g.eKE/g.nSph));
   }
}
#endif

inline double pressSpec(Cparticle &p,CglobalVars &g) {
   if (PERIODIC[0]) {
      return (p.press);
   } else {
      return (p.press-(g.custom[8]/4));
   }
}

inline double viscSpecAbs(Cparticle &p,CglobalVars &g) {
   if (PERIODIC[0]) {
      return dot(p.fv,p.fv);
   } else {
      return (dot(p.fv,p.fv)-g.custom[7]);
   }
}

inline double viscSpec(Cparticle &p,CglobalVars &g) {
   if (PERIODIC[0]) {
      return dot(p.fv,p.vhat);
   } else {
      return (dot(p.fv,p.vhat)-g.custom[1]);
   }
}           

inline double enstrophySpec(Cparticle &p,CglobalVars &g) {
   if (PERIODIC[0]) {
      return len2(p.vort);
   } else {
      return (len2(p.vort)-2.0*g.custom[0]);
   }
}

inline void zeroV(Cparticle &p,CglobalVars &g) {
   p.v = 0;
   p.vhat = 0;
   g.eKE = 0;
}

inline void setVelSpec(Cparticle &p,CglobalVars &g,vectInt &coords) {
   const double mult = 2.0*PI/(RMAX[0]-RMIN[0]);
   const double k2 = len2(coords);
   double c;
   if (k2>0) {
      if (k2<8*8) {
         c = 6.25*pow(k2,-5.0/6.0);
      } else {
         c = 100*pow(k2,-3.0/2.0);
      }
      //c = 0.1;
      p.v[0] += c*sin(mult*coords[0]*(p.r[0]-RMIN[0])+ mult*coords[1]*(p.r[1]-RMIN[1]));
      p.vhat = p.v;
      g.eKE += 0.5*p.mass*len2(p.v);
   }
}

inline double velxSpec(Cparticle &p,CglobalVars &g) {
   return p.v[0];
}

void CcustomOutput::calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io, CsphIncompress *sph) {

   //char strTimestep[20];
   //sprintf(strTimestep,"%7.7d",data->globals.outstep);
   //string outputFN = io->getFilename()+filename+strTimestep+".dat";
   //string filename = "maxV";
   //string outputFN = io->getFilename()+filename+".dat";
   //ofstream newFile(outputFN.c_str(),ios::out|ios::app);
   //newFile <<"#initial_dist velStruct"<<endl;
   //data->globals.maxV = 0;
   //data->globals.aveDens = 0;
   //data->globals.varDens = 0;
   //data->traverse<calcMaxVandDens,ifSphOrSphBoundary>();
   //data->globals.n = data->getParticles()->size();
   //double tmpData[2];
   //tmpData[0] = data->globals.n;
   //tmpData[1] = data->globals.aveDens;
   //data->sumOverProcs(tmpData,2);
   //data->globals.aveDens = tmpData[1]/tmpData[0];
   //data->traverse<calcTheVarDens,ifSphOrSphBoundary>();
   //data->sumOverProcs(&(data->globals.varDens),1);
   //data->globals.varDens /= tmpData[0];
   //newFile<<data->globals.sphStep<<' '<<data->globals.outstep<<' '<<data->globals.time<<' '<<data->globals.maxV<<' '<<data->globals.aveDens<<' '<<data->globals.varDens<<endl;

   //data->traverse<zeroV,ifSph>();
   //const int numModes = int(0.4*(RMAX[1]-RMIN[1])/PSEP);
   //setSpectrum<ifSph,setVelSpec>(data,numModes);

   //if (data->globals.mpiRank==0) cout<<"\tRe-initialising density..."<<endl;
   //for (int i=0;i<10;i++) {
   //   data->reset();
   //   data->traverse<calcHFromDens,ifSphOrSphBoundary>();
   //   data->neighboursGroup<calcDensity,ifSphOrSphBoundary>();
   //}


   //data->reset();

/*
   vector<Cparticle> psGrid;
   vectInt gridDims;
   sph->renderToGrid(psGrid,gridDims);
   
   ofstream fo;
   string outputFN = "VGrid";
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   outputFN = io->getFilename()+outputFN+strTimestep+".dat";

   cout <<"writing grid data to file: "<<outputFN<<endl;
   fo.open(outputFN.c_str(),ios::out|ios::trunc);
   fo << "# "<<NX<<" x "<<NY<<" matrix"<<endl;

   for (particleContainer::iterator p = psGrid.begin();p != psGrid.end();p++) {
      if (p->r[0] == RMIN[0]) {
         fo<<endl;
      }
      fo<<len(p->v)<<' ';
   }
   fo.close();
*/

   //io->writeGrid(data->globals.outstep,psGrid,gridDims,&(data->globals));
          
   //calcVelocityStructure<ifInteriorSph>(data,rng,io);
   //calcFullVelocityStructure<ifInteriorSph>(data,io);
   
   //calcPositionStructure<ifInteriorSph>(data,rng,io);

   //pairDoubling(data,pairData,openFiles,rng,io);
   
   //outputVelocitySpectrum<ifSph>(data,io);

   //string name; 
   //name = "VelSpec";
   //outputSpectrum<ifSph,velxSpec>(data,io,name);
  
   string name; 
   //name = "PressSpec";
   //outputSpectrum<ifSph,pressSpec>(data,io,name);
#ifdef AVE_VELOCITY
   name = "aveVSpec";
   outputSpectrum<ifSph,aveVSpec>(data,io,name);
#endif
   name = "VelSpec";
   outputSpectrum<ifSph,velSpec>(data,io,name);
   //name = "ViscSpec";
   //outputSpectrum<ifSph,viscSpec>(data,io,name);
   //name = "ViscSpecAbs";
   //outputSpectrum<ifSph,viscSpecAbs>(data,io,name);
   //name = "EnstrophySpec";
   //outputSpectrum<ifSph,enstrophySpec>(data,io,name);

}
