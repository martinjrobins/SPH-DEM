#include "customOutput.h"
#include "sph.h"
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
   io->writeGrid(outstep,ps,gridDims);
#endif
#ifdef PAIR_DBL
   //pair doubling stuff
   //go through existing pairs
   Sgi::slist<CpairData>::iterator i = pairData.begin();
   Sgi::slist<CpairData>::iterator i_before;
   while (i != pairData.end()) {
      double dt = data->globals.time-i->initT;
      double r2 = len2(i->p1->r-i->p2->r);
      if ((r2>i->finalR2)||(dt>PAIR_DBL_MAX_TIME)) {
         double lyap;
         if (dt>PAIR_DBL_MAX_TIME) {
            lyap = 0;
         } else {
            lyap = (0.5/dt) * log(r2/i->initR2);
         }
         *(i->outFile) << sqrt(i->initR2) <<" "<< lyap <<" "<<sqrt(r2)<<" "<<dt<<endl;
         if (i==pairData.begin()) {
            i_before = i;
            i++;
            pairData.erase(i_before);
         } else {
            i++;
            pairData.erase_after(i_before);
         }
      } else {
         i_before = i;
         i++;
      }
   }

   //open new output file
   string outputFN = "PairDoublingTimes";
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",outstep);
   outputFN = io->getFilename()+outputFN+strTimestep+".dat";
   ofstream *newFile = new ofstream(outputFN.c_str());
   *newFile <<"#initial_dist lyaponov_exponent final_dist time_taken"<<endl;
   openFiles.push_back(newFile);
   
   //init new pairs
   vector<Cparticle> *ips = data->getParticles();
   int numParticles = ips->size();
   for (int i=0;i<PAIR_DBL_SAMPLES_PER_OUTSTEP;i++) {
      int ip1 = (int)floor(gsl_rng_uniform(rng)*numParticles);
      while ((*ips)[ip1].iam != sph) {
         ip1 = (int)floor(gsl_rng_uniform(rng)*numParticles);
      }
      int bin = (int)ceil(gsl_rng_uniform(rng)*PAIR_DBL_BINS);
      double r2 = pow(0.5*(RMAX[0]-RMIN[0])/pow(2.0,PAIR_DBL_BINS-bin),2);
      double r2_floor = pow(0.5*(RMAX[0]-RMIN[0])/pow(2.0,PAIR_DBL_BINS-bin+0.5),2);
      double r2_ceil = pow(0.5*(RMAX[0]-RMIN[0])/pow(2.0,PAIR_DBL_BINS-bin-0.5),2);
      vector<Cparticle *> candidates;
      for (int j=0;j<numParticles;j++) {
         double thisR2 = len2((*ips)[j].r-(*ips)[ip1].r);
         if ((thisR2>=r2_floor)&&(thisR2<=r2_ceil)&&((*ips)[j].iam==sph)) {
            candidates.push_back(&(*ips)[j]);
         }
      }
      
      if (!candidates.empty()) {
         int ip2 = (int)floor(gsl_rng_uniform(rng)*candidates.size());
         CpairData newPD;
         newPD.p1 = &((*ips)[ip1]);
         newPD.p2 = candidates[ip2];
         newPD.outFile = newFile;
         newPD.initT = data->globals.time;
         newPD.initR2 = len2(newPD.p1->r-newPD.p2->r);
         newPD.finalR2 = newPD.initR2*PAIR_DBL_R_MULT2;
         pairData.push_front(newPD);
      }
   }

   cout <<"\t\tThere are now "<<pairData.size()<<" pairs being tracked."<<endl;
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
