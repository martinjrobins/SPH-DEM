#include "customOutput.h"
#include "sph.h"
#include "misc.impl.h"
#include "dataLL.impl.h"

using namespace Nsph;

inline void calcEnstrophyAndAveDens(Cparticle &p,CglobalVars &g) {
   g.custom[0] += 0.5*p.mass*pow(p.vort,2)/p.dens;
   g.custom[1] += p.dens;
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

#ifdef FFT
   //setup a grid of particles over the box
   vector<Cparticle> ps;
   ps.resize((NX+1)*(NY+1));
   double cang = cos(custSim->getFrameAng());
   double sang = sin(custSim->getFrameAng());
   for (int i=0;i<=NX;i++) {
      for (int j=0;j<=NY;j++) {
         int ind = i*(NY+1)+j;
         ps[ind].tag = ind;
         vect rtmp;
         rtmp[0] = i*PSEP+RMIN[0];
         rtmp[1] = j*PSEP+RMIN[1];
         ps[ind].r[0] = rtmp[0]*cang-rtmp[1]*sang;
         ps[ind].r[1] = rtmp[0]*sang+rtmp[1]*cang;
      }
   }
   vectInt gridDims;
   gridDims[0] = NX+1;
   gridDims[1] = NY+1;
   data->neighboursUsing<calcVonGrid>(ps);
   // perform fft TODO
    
   // fftw_execute(p);
   io->writeGrid(outstep,ps,gridDims);
#endif

#ifdef VORT_LEASTSQUARES
   data->neighboursGroup<calcVortLeastSquares,ifSph>();
#else
   data->neighboursGroup<calcVortSPHSum,ifSph>();
#endif
   data->globals.custom[0] = 0;
   data->globals.custom[1] = 0;
   data->traverse<calcEnstrophyAndAveDens,ifSphOrSphBoundary>();
   data->globals.custom[1] /= data->numParticles();
}

inline bool ifInteriorSph(Cparticle &p) {
   return (p.iam == sph)&&(abs(p.r[0])<RMAX[0]*0.75)&&(abs(p.r[1])<RMAX[1]*0.75);
}

void CcustomOutput::calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io) {

   string velStructFilename = "velStruct";
   calcVelocityStructure<ifInteriorSph>(data,velStructFilename,rng,io);
   
   string posStructFilename = "posStruct";
   calcPositionStructure<ifInteriorSph>(data,posStructFilename,rng,io);
 

   //pairDoubling(data,pairData,openFiles,rng,io);
   /*
   vector<double> *oneD = calcVelocitySpectrum<ifInteriorSph>(data);
   ofstream fo;
   string outputFN = "VelSpec";
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   outputFN = io->getFilename()+outputFN+strTimestep+".dat";
   cout <<"writing 1D spectrum to file: "<<outputFN<<endl;
   fo.open(outputFN.c_str());
   fo << "# mode 1DSpectrum"<<endl;
   int n = oneD->size();
   for (int i=0;i<n;i++) {
      fo<<i<<' '<<(*oneD)[i]<<endl;
   }
   fo.close();
   delete oneD;
   */
   cout <<"finished post process"<<endl;
}
