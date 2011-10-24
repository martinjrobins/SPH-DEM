
#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H
#include "vect.h"
#include <vector>
#include <list>
#include "sphConstants.h"
#include "demConstants.h"
#include <blitz/array.h>

class sphGhost:baseGhost {
public:
#ifdef SMOOTHED_VISC_VELOCITY
      vect viscv;
#endif
#ifdef AVE_VELOCITY
      vect aveV;
#endif
#ifdef SLK
      vect currR;
      vect dr[SLK_NUM_TIMESTEPS];
      //double drho[SLK_NUM_TIMESTEPS];
      vect thisRInc;
#endif
#ifdef LIQ_DEM
      double porosity;
#endif
      bool concave;
      vect norm1;
      vect norm2;
#ifdef _3D_
      vect norm3;
#endif
      vect vort;
      sphGhost& operator=(const sphParticle &p);
};

class sphParticle:baseParticle {
   public:
   typedef sphGhost ghost; 
   sphParticle() {
      baseParticle();
#ifdef CORRECTED_GRADIENT
      invM = 0;
      sum = 0;
      sumGrad = 0;
#endif
#ifdef LIQ_DEM
      porosity = 1.0;
      fdrag = 0.0;
#endif
      vort = 0;
#ifdef SMOOTHED_VISC_VELOCITY
      viscv = 0;
#endif
#ifdef AVE_VELOCITY
      aveV = 0;
      for (int i=0;i<AVE_VELOCITY_N;i++) {
         pastVs[i] = 0;
      }
#endif
#ifdef SLK
      for(int i=0;i<SLK_NUM_TIMESTEPS;i++) {
         dr[i] = 0;
         drho[i] = 0;
      }
      currR = r;
      thisRInc = 0;
      dRhoKernel = 0;
      dMassDiff = 0;
#endif
      vhat0 = 0;
      v0 = 0;
#ifdef SAVE_VHALF
      vHalf = 0;
      rHalf = 0;
#endif
      fp = 0;
      fv = 0;
      fg = 0;
      dddt = 0;
      press = 0;
      pdr2 = 0;
      spsound = 0;
      alpha = 1;
      u = 0;
      vort = 0;
      dudt = 0;
      colour = 0;
      eViscF = 0;
      eViscB = 0;
      eBForce = 0;
      eFF = 0;
      norm1 = 0;
      norm2 = 0;
#ifdef _3D_
      norm3 = 0;
#endif
      concave = true;
      dist = 0;

#ifdef SPH_SMOOTHING
      matA.resize(NUM_SMOOTH_NEIGHBOURS);
#endif
   };

   sphParticle(const sphParticle &p) {
      baseParticle(p);
#ifdef SPH_SMOOTHING
      matA.resize(NUM_SMOOTH_NEIGHBOURS);
#endif
   }
      
   sphParticle& operator=(const sphParticle &p) {
      static_cast<baseParticle>(*this) = p;
#ifdef LIQ_DEM
      porosity = p.porosity;
#endif
      v0 = p.v0;
      vort = p.vort;
#ifdef SAVE_VHALF
      vHalf = p.vHalf;
      rHalf = p.rHalf;
#endif
      vhat0 = p.vhat0;
#ifdef SMOOTHED_VISC_VELOCITY
      viscv = p.viscv;
#endif
#ifdef AVE_VELOCITY
      aveV = p.aveV;
      for (int i=0;i<AVE_VELOCITY_N;i++) {
         pastVs[i] = p.pastVs[i];
      }
#endif
#ifdef SLK
      for(int i=0;i<SLK_NUM_TIMESTEPS;i++) {
         dr[i] = p.dr[i];
         drho[i] = p.drho[i];
      }
      currR = p.currR;
      thisRInc = p.thisRInc;
      dRhoKernel = p.dRhoKernel;
      dMassDiff = p.dMassDiff;
#endif
      fp = p.fp;
      fv = p.fv;
      fg = p.fg;
      fb = p.fb;
      ff = p.ff;
      norm1 = p.norm1;
      norm2 = p.norm2;
#ifdef _3D_
      norm3 = p.norm3;
#endif
      concave = p.concave;
      dist = p.dist;
      vort = p.vort;
      dddt = p.dddt;
      press = p.press;
      pdr2 = p.pdr2;
      spsound = p.spsound;
      alpha = p.alpha;
      dudt = p.dudt;
      u = p.u;
      colour = p.colour;

      eViscF = p.eViscF;
      deViscFdt = p.deViscFdt;
      eViscB = p.eViscB;
      deViscBdt = p.deViscBdt;
      eBForce = p.eBForce;
      deBForcedt = p.deBForcedt;
      eFF = p.eFF;
      deFFdt = p.deFFdt;
      return *this;
   }
   
#ifdef CORRECTED_GRADIENT
   vectTensor invM;
   double sum;
   vect sumGrad;
#endif
#ifdef LIQ_DEM
   double porosity;
   vect dVortDt;
   vect eta[3];
   vect tdv[3];
   int gHead;
   vect fdrag;
   vect gVect[TWIN_N];
#endif
   vect v0;
#ifdef SAVE_VHALF
   vect vHalf;
   vect rHalf;
#endif
   vect vhat0;
   vect fp,fv,fg,fb,ff;

   vect norm1,norm2;
#ifdef _3D_
   vect norm3;
#endif
   double dist;
   bool concave;
   vect vort;
   double dddt;
   double press;
   double pdr2;
   double spsound;
   double alpha;
   double dudt;
   double u;
   double colour;
#ifdef SMOOTHED_VISC_VELOCITY
   vect viscv;
#endif
#ifdef SLK
   vect dr[SLK_NUM_TIMESTEPS];
   double drho[SLK_NUM_TIMESTEPS];
   vectTensor G;
   vect currR;
   vect thisRInc;
   double dRhoKernel;
   double dMassDiff;
#endif
#ifdef AVE_VELOCITY
   vect aveV;
   vect pastVs[AVE_VELOCITY_N];
#endif
#ifdef DIRECT_SMOOTHING
   vect fhat;
#endif
#ifdef SPH_SMOOTHING
   vect residual,basis;
   vect pTimesA;
   vector<double> matA;
#endif
   double eViscF,deViscFdt;
   double eViscB,deViscBdt;
   double eBForce,deBForcedt;
   double eFF,deFFdt;
   double tmp;

   #ifdef LIQ_DEM
   RETVARIABLE(porosity,double,sphParticle,p.porosity)
   RETVARIABLE(fdrag,double,sphParticle,p.fdrag)
   #endif
   RETVARIABLE(fp,double,sphParticle,p.fp)
   RETVARIABLE(fv,double,sphParticle,p.fv)
   RETVARIABLE(fb,double,sphParticle,p.fb)
   RETVARIABLE(ff,double,sphParticle,p.ff)
   RETVARIABLE(vorticity,vect,sphParticle,p.vort)
   RETVARIABLE(colour,double,sphParticle,p.colour)
};

inline sphGhost& sphGhost::operator=(const sphParticle &p) {
      static_cast<baseGhost>(*this) = p;
#ifdef LIQ_DEM
         porosity = p.porosity;
#endif
         vort = p.vort;
#ifdef SMOOTHED_VISC_VELOCITY
         viscv = p.viscv;
#endif
#ifdef AVE_VELOCITY
         aveV = p.aveV;
#endif
#ifdef SLK
         for(int i=0;i<SLK_NUM_TIMESTEPS;i++) {
            dr[i] = p.dr[i];
            //drho[i] = p.drho[i];
         }
         currR = p.currR;
         thisRInc = p.thisRInc;
#endif
         norm1 = p.norm1*(p.dist+1);
         norm2 = p.norm2;
#ifdef _3D_
         norm3 = p.norm3;
#endif
         concave = p.concave;
         return *this;
};

inline sphParticle& sphParticle::operator=(const sphGhost &g) {
      static_cast<baseParticle>(*this) = g;
#ifdef LIQ_DEM
      porosity = g.porosity;
#endif
      vort = g.vort;
#ifdef SMOOTHED_VISC_VELOCITY
      viscv = g.viscv;
#endif
#ifdef AVE_VELOCITY
      aveV = g.aveV;
#endif
#ifdef SLK
      for(int i=0;i<SLK_NUM_TIMESTEPS;i++) {
        dr[i] = g.dr[i];
        //drho[i] = g.drho[i];
      }
      currR = g.currR;
      thisRInc = g.thisRInc;
#endif
      double r1 = len(g.norm1);
      norm1 = g.norm1/r1;
      norm2 = g.norm2;
#ifdef _3D_
      norm3 = g.norm3;
#endif
      dist = r1-1;
      concave = g.concave;
      return *this;
};


