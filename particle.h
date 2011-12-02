
#ifndef PARTICLE_H
#define PARTICLE_H
#include "vect.h"
#include <vector>
#include <list>
#include "sphConstants.h"
#include "demConstants.h"
#include <blitz/array.h>

typedef vect my2dMatrix[NDIM]; 

/*
 * enum type of different particles types
 * used for the iam member of the Cparticle type
 */
enum iamTypes {sph, sphBoundary, immovable, boundary,ghost,boundaryobj1,boundaryobj2,boundaryobj3,boundaryobj4,boundaryobj5,boundaryobj6,dem,demBoundary};


class Cparticle;
typedef vector<Cparticle> particleContainer;

/*
 * class type for ghost data (the data that is passed between cpus instead of
 * the full particle data)
 */
class CghostData {
public:
      vect r;
      vect v;
      vect vhat;
#ifdef SLK
      vect currR;
      vect dr[SLK_NUM_TIMESTEPS];
      vect thisRInc;
#endif
#ifdef LIQ_DEM
      double porosity;
#endif
      double mass;
      double dens;
      double h;
      bool concave;
      vect norm1;
      vect norm2;
#ifdef _3D_
      vect norm3;
#endif
      vect vort;
      int tag;
      enum iamTypes iam;

      CghostData& operator=(const Cparticle &p);
};


/*
 * main particle class. Contains all useful data on each particle
 */
class Cparticle {
   public:

   /*
    * long constructor. Contains all the default values for all of the
    * variables.
    */
   Cparticle() {
#ifdef CORRECTED_GRADIENT
      invM = 0;
      sumGrad = 0;
      sum = 0;
#endif
#ifdef VAR_H_CORRECTION
      omega = 0;
#endif
#ifdef LIQ_DEM
      porosity = 1.0;
      dporositydt = 0.0;
      gHead = 0;
      fdrag = 0.0;
      for (int i=0;i<TWIN_N;i++) {
         gVect[i] = 0.0;
      }
      fvel = 0.0;
      shepSum = 0.0;
#endif
      r = 0;
      v = 0;
      vhat = 0;
      vort = 0;
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
      fp = 0;
      fv = 0;
      fg = 0;
      f = 0;
      h = 0;
      dens = 0;
      tmp = 0;
      dddt = 0;
      press = 0;
      pdr2 = 0;
      u = 0;
      vort = 0;
      dudt = 0;
      iam = sph;
      tag = 0;
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
   };

   /*
    * init particle with starting position
    */
   Cparticle(vect init) {
      Cparticle();
      r = init;
   };

   /*
    * copy constructor
    */
   Cparticle(const Cparticle &p) {
      *this = p;
   }
      
   /*
    * copy a particle onto another. 
    */
   Cparticle& operator=(const Cparticle &p) {
#ifdef LIQ_DEM
      porosity = p.porosity;
      dporositydt = p.dporositydt;
      gHead = p.gHead;
      fdrag = p.fdrag;
      for (int i=0;i<TWIN_N;i++) {
         gVect[i] = p.gVect[i];
      }
      fvel = p.fvel;
      shepSum = p.shepSum;
#endif
      mass = p.mass;
      r = p.r;
      v = p.v;
      v0 = p.v0;
      vort = p.vort;
      vhat = p.vhat;
      vhat0 = p.vhat0;
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
      f = p.f;
      fb = p.fb;
      ff = p.ff;
      norm1 = p.norm1;
      norm2 = p.norm2;
#ifdef _3D_
      norm3 = p.norm3;
#endif
      concave = p.concave;
      vort = p.vort;
      h = p.h;
      dens = p.dens;
      tmp = p.tmp;
      dddt = p.dddt;
      press = p.press;
      pdr2 = p.pdr2;
      minDt = p.minDt;
      dudt = p.dudt;
      u = p.u;
      iam = p.iam;
      tag = p.tag;
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
   

      
   Cparticle& operator=(const CghostData &g);

   /*
    * definition of a neighbour
    */
   bool is_neighbr(Cparticle &_p) {
      return (&_p != this) && (len2(r-_p.r) <= pow(KERNAL_RADIUS*max(h,_p.h),2));
   }

   /*
    * list of particle variables
    */
#ifdef CORRECTED_GRADIENT
   vectTensor invM;
   vect sumGrad;
   double sum;
#endif
#ifdef VAR_H_CORRECTION
   double omega;
#endif
#ifdef LIQ_DEM
   double porosity;
   double dporositydt;
   int gHead;
   vect fdrag;
   vect gVect[TWIN_N];
   vect fvel;
   double shepSum;
#endif
   double mass;
   vect r;
   vect v;
   vect v0;
   vect vhat,vhat0;
   vect fp,fv,fg,f,fb,ff;

   vect norm1,norm2;
#ifdef _3D_
   vect norm3;
#endif
   bool concave;
   vect vort;
   double h;
   double dens,dddt;
   double tmp;
   double press;
   double pdr2;
   double minDt;
   double dudt;
   double u;
   enum iamTypes iam;
   int tag;
   double colour;
#ifdef SLK
   vect dr[SLK_NUM_TIMESTEPS];
   double drho[SLK_NUM_TIMESTEPS];
   vectTensor G;
   vect currR;
   vect thisRInc;
   double dRhoKernel;
   double dMassDiff;
#endif

   double eViscF,deViscFdt;
   double eViscB,deViscBdt;
   double eBForce,deBForcedt;
   double eFF,deFFdt;

};

/*
 * these assignment operators are important as they map particles onto ghost
 * particles and visa versa.
 */
inline CghostData& CghostData::operator=(const Cparticle &p) {
#ifdef LIQ_DEM
         porosity = p.porosity;
#endif
         r = p.r;
         v = p.v;
         vort = p.vort;
         vhat = p.vhat;
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
         mass = p.mass;
         dens = p.dens;
         h = p.h;
         tag = p.tag;
         iam = p.iam;
         norm1 = p.norm1;
         norm2 = p.norm2;
#ifdef _3D_
         norm3 = p.norm3;
#endif
         concave = p.concave;
         return *this;
};

inline Cparticle& Cparticle::operator=(const CghostData &g) {
#ifdef LIQ_DEM
      porosity = g.porosity;
#endif
      r = g.r;
      v = g.v;
      vhat = g.vhat;
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
      }
      currR = g.currR;
      thisRInc = g.thisRInc;
#endif
      h = g.h;
      mass = g.mass;
      dens = g.dens;
      tag = g.tag;
      double r1 = len(g.norm1);
      norm1 = g.norm1;
      norm2 = g.norm2;
#ifdef _3D_
      norm3 = g.norm3;
#endif
      concave = g.concave;
      iam = g.iam;
      return *this;
};


#endif
