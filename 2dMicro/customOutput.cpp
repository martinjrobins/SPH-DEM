#include "customOutput.h"
#include "sph.h"
#include "misc.impl.h"
#include "dataLL.impl.h"

using namespace Nsph;

inline void calcEnstrophyAndAveDens(Cparticle &p,CglobalVars &g) {
   g.custom[0] += 0.5*p.mass*pow(p.vort,2)/p.dens;
   g.custom[1] += p.mass*dot(p.vhat,p.fv)/p.dens;
   g.custom[2] += p.mass*dot(p.vhat,p.ff)/p.dens;
   g.custom[3] += p.mass*dot(p.vhat,p.f)/p.dens;
   g.custom[4] += p.mass*p.pdr2*p.dddt/p.dens;
   g.custom[5] += p.mass*dot(p.vhat,p.fp)/p.dens;
   g.custom[6] += 0.5*g.dt*p.mass*dot(p.ff,p.ff)/p.dens;
   g.custom[7] += p.mass*dot(p.fv,p.fv)/p.dens;
   //g.custom[8] += p.mass*dot(p.f,p.vhat)/p.dens;
}

inline void calcVonGrid(Cparticle &p,vector<Cparticle *> &neighbrs,CglobalVars &g) {
   double b0;
   vect bRest;
   vector<double> vWab;
   calcB_MLS(p,neighbrs,b0,bRest,vWab);
   int n = neighbrs.size();
   p.v = 0;
   p.vort = 0;
   for (int i=0;i<n;i++) {
      Cparticle *pn = neighbrs[i];
      p.v += (pn->v*pn->mass/pn->dens)*W_MLS(p.r-pn->r,vWab[i],b0,bRest);
      p.vort += (pn->vort*pn->mass/pn->dens)*W_MLS(p.r-pn->r,vWab[i],b0,bRest);
   }
}


void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {
/*
   //setup a grid of particles over the box
   vector<Cparticle> ps;
   ps.resize((NX+1)*(NY+1));
   for (int i=0;i<=NX;i++) {
      for (int j=0;j<=NY;j++) {
         int ind = i*(NY+1)+j;
         ps[ind].tag = ind;
         vect rtmp;
         rtmp[0] = i*PSEP+RMIN[0];
         rtmp[1] = j*PSEP+RMIN[1];
         ps[ind].r[0] = rtmp[0];
         ps[ind].r[1] = rtmp[1];
      }
   }
   vectInt gridDims;
   gridDims[0] = NX+1;
   gridDims[1] = NY+1;
   data->neighboursUsing<calcVonGrid>(ps);
   io->writeGrid(outstep,ps,gridDims);
#ifdef VORT_LEASTSQUARES
   data->neighboursGroup<calcVortLeastSquares,ifSph>();
#else
   data->neighboursGroup<calcVortSPHSum,ifSph>();
#endif
*/
   data->globals.custom[0] = 0;
   data->globals.custom[1] = 0;
   data->globals.custom[2] = 0;
   data->globals.custom[3] = 0;
   data->globals.custom[4] = 0;
   data->globals.custom[5] = 0;
   data->globals.custom[6] = 0;
   data->globals.custom[7] = 0;
   //data->globals.custom[8] = 0;
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


inline void zeroV(Cparticle &p,CglobalVars &g) {
   p.v = 0;
   p.vhat = 0;
   g.eKE = 0;
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
  
   /*string name; 
   name = "VelSpec";
   outputSpectrum<ifSph,velSpec>(data,io,name);
   name = "ViscSpec";
   outputSpectrum<ifSph,viscSpec>(data,io,name);
   name = "ViscSpecAbs";
   outputSpectrum<ifSph,viscSpecAbs>(data,io,name);
   name = "EnstrophySpec";
   outputSpectrum<ifSph,enstrophySpec>(data,io,name);
   */
}
