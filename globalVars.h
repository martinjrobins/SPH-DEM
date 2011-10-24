#ifndef GLOBALVARS_H
#define GLOBALVARS_H

#include "vect.h"
#include <vector>
#include <blitz/array.h>
#include <sys/time.h>

#define GLOBAL_CUSTOM_BUFFER_SIZE 9

class CglobalVars {
   public:
      CglobalVars() {
         sphStep = 0;
         time = 0;
         dt = 0;
         newDt = 0;
         angMom = 0;
         linMom = 0;
         eElast = 0;
         eElastExact = 0;
         eKE = 0;
         eFF = 0;
         eViscF = 0;
         eViscB = 0;
         eBForce = 0;
         eTotal = 0;

         edElastdt = 0;
         edFFdt = 0;
         edViscFdt = 0;
#ifdef SLK
         slkIndex0 = 0;
         slkIndexLast = SLK_NUM_TIMESTEPS-1;
#endif

         maxdt = 0;
         procNeighbrs.resize(3);
         procNeighbrs = -1;
         for (int i=0;i<NDIM;i++) {
            procDomain[2*i] = RMIN[i];
            procDomain[2*i+1] = RMAX[i];
         }

         mpiSize = 1;
         mpiRank = 0;
         n = 0;
         nSph = 0;
         maxV = 0;
         maxF = 0;
         maxFF = 0;
      };

      CglobalVars(const CglobalVars &g) {
         angMom = g.angMom;
         linMom = g.linMom;
         eBForce = g.eBForce;
         eElast = g.eElast;
         eElastExact = g.eElastExact;
         eKE = g.eKE;
         eViscF = g.eViscF;
         eViscB = g.eViscB;
         eTotal = g.eTotal;
         eFF = g.eFF;

         edElastdt = g.edElastdt;
         edFFdt = g.edFFdt;
         edViscFdt = g.edViscFdt;

         dt = g.dt;
         newDt = g.newDt;
         time = g.time;
         sphStep = g.sphStep;
         for (int i=0;i<GLOBAL_CUSTOM_BUFFER_SIZE;i++) {
            custom[i] = g.custom[i];
         }
         mpiRank = g.mpiRank;
         mpiSize = g.mpiSize;
         n = g.n;
         nSph = g.nSph;
         maxV = g.maxV;
         maxF = g.maxF;
         maxFF = g.maxFF;
         aveDens = g.aveDens;
         aveDensFromMass = g.aveDensFromMass;
         varDens = g.varDens;
         rmsFF = g.rmsFF;

         wtTotalMPI = g.wtTotalMPI;
         wtTotalUpdateDomain = g.wtTotalUpdateDomain;
         wtTotalOutStep = g.wtTotalOutStep;
         wtTotalForceCalc = g.wtTotalForceCalc;
         wtTotalDddtCalc = g.wtTotalDddtCalc;
         wtTotalReset = g.wtTotalReset;
#ifdef SMOOTHING
         wtTotalSmooth = g.wtTotalSmooth;
         wtTotalSmoothMPI = g.wtTotalSmoothMPI;
#endif
      }


      CglobalVars& operator+=(CglobalVars &g) {
         angMom += g.angMom;
         linMom += g.linMom;
         eBForce += g.eBForce;
         eElast += g.eElast;
         eElastExact += g.eElastExact;
         eKE += g.eKE;
         eViscF += g.eViscF;
         eViscB += g.eViscB;
         eFF += g.eFF;
         eTotal += g.eTotal;

         edElastdt += g.edElastdt;
         edFFdt += g.edFFdt;
         edViscFdt += g.edViscFdt;

         maxV = max(maxV,g.maxV);
         maxF = max(maxF,g.maxF);
         maxFF = max(maxFF,g.maxFF);

         for (int i=0;i<GLOBAL_CUSTOM_BUFFER_SIZE;i++) {
            custom[i] += g.custom[i];
         }
         return *this;
      }

      void init(CglobalVars &g) {
         time = g.time;
         sphStep = g.sphStep; 
         dt = g.dt;
         newDt = g.newDt;
         maxdt = g.maxdt;
         mpiSize = g.mpiSize;
         aveDens = g.aveDens;
         aveDensFromMass = g.aveDensFromMass;
         varDens = g.varDens;
         rmsFF = g.rmsFF;
         n = g.n;
         nSph = g.nSph;

         wtTotalMPI = g.wtTotalMPI;
         wtTotalUpdateDomain = g.wtTotalUpdateDomain;
         wtTotalOutStep = g.wtTotalOutStep;
         wtTotalForceCalc = g.wtTotalForceCalc;
         wtTotalDddtCalc = g.wtTotalDddtCalc;
         wtTotalReset = g.wtTotalReset;
#ifdef SMOOTHING
         wtTotalSmooth = g.wtTotalSmooth;
         wtTotalSmoothMPI = g.wtTotalSmoothMPI;
#endif
      }

      void startOutStep() {
         wtTotalMPI.tv_sec = 0;
         wtTotalMPI.tv_usec = 0;
         wtTotalUpdateDomain.tv_sec = 0;
         wtTotalUpdateDomain.tv_usec = 0;
         gettimeofday(&wtTotalOutStep,NULL);
         wtTotalForceCalc.tv_sec = 0;
         wtTotalForceCalc.tv_usec = 0;
         wtTotalDddtCalc.tv_sec = 0;
         wtTotalDddtCalc.tv_usec = 0;
         wtTotalReset.tv_sec = 0;
         wtTotalReset.tv_usec = 0;
#ifdef SMOOTHING
         wtTotalSmooth.tv_sec = 0;
         wtTotalSmooth.tv_usec = 0;
         wtTotalSmoothMPI.tv_sec = 0;
         wtTotalSmoothMPI.tv_usec = 0;
#endif
      }

      void endOutStep() {
         timeval tmp;
         gettimeofday(&tmp,NULL);
         wtTotalOutStep.tv_sec = tmp.tv_sec-wtTotalOutStep.tv_sec;
         wtTotalOutStep.tv_usec = tmp.tv_usec-wtTotalOutStep.tv_usec;
      }

      double angMom;
      vect linMom;
      double eBForce,eElast,eElastExact,eKE,eViscF,eViscB,eFF,eTotal;
      double edElastdt,edViscFdt,edFFdt;
      double maxV, aveDens, varDens, aveDensFromMass;
      double maxF, maxFF,rmsFF;
      double dt,newDt,maxdt;
      int n,nSph;
#ifdef SPH_SMOOTHING
      vect alpha,beta,residual2,alphaDenom;
#endif
#ifdef DIRECT_SMOOTHING
      vect residual;
      int DSTmpInt;
#endif
#ifdef SMOOTHING
      timeval wtTotalSmooth,wtTotalSmoothMPI;
#endif
#ifdef GRID_SMOOTHING
      double itError2;
#endif
#ifdef SLK
      unsigned int slkIndex0;
      unsigned int slkIndexLast;
#endif

      double time;
      int sphStep;
      double custom[GLOBAL_CUSTOM_BUFFER_SIZE];
      double procDomain[NDIM*2];
      Array<int,NDIM> procNeighbrs;
      int mpiRank;
      int mpiSize;

      double tmp;
      double tmpArray[5];

      int outstep;

      
      timeval wtTotalMPI,wtTotalUpdateDomain,wtTotalOutStep,wtTotalForceCalc,wtTotalDddtCalc,wtTotalReset;
    
};

#endif
