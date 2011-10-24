
#ifndef BASEPARTICLE_H
#define BASEPARTICLE_H
#include "vect.h"
#include <vector>
#include <list>
#include "sphConstants.h"
#include "demConstants.h"
#include <blitz/array.h>

#define RETVARIABLE(name, type, particleType, funct)                            \
   class name {                                                                 \
      public:                                                                   \
        type& operator()(particleType p) {                                      \
           return funct;                                                        \
        }                                                                       \
        const string = 'name';                                                  \
   }                                                                            \



#ifdef SPH_SMOOTHING
const int NUM_SMOOTH_NEIGHBOURS = int(pow(HFAC*SMOOTH_HFAC_MULT*4.0+1,NDIM));
#endif

typedef vect my2dMatrix[NDIM]; 

enum iamTypes {sph, sphBoundary, immovable, boundary,ghost,boundaryobj1,boundaryobj2,boundaryobj3,boundaryobj4,boundaryobj5,boundaryobj6,dem};


class baseGhost {
public:
      vect r;
      vect v;
      vect vhat;
      double mass;
      double dens;
      double h;
      int tag;
      enum iamTypes iam;
      CghostData& operator=(const Cparticle &p);
};


class baseParticle {
   public:
   typedef baseGhost ghost; 
   baseParticle() {
      r = 0;
      v = 0;
      vhat = 0;
      f = 0;
      h = 0;
      dens = 0;
      iam = sph;
      tag = 0;
   };
   baseParticle(vect init) {
      baseParticle();
      r = init;
   };

   baseParticle(const baseParticle &p) {
      *this = p;
   }

   baseParticle& operator=(const baseParticle &p) {
      mass = p.mass;
      r = p.r;
      v = p.v;
      vhat = p.vhat;
      f = p.f;
      h = p.h;
      dens = p.dens;
      minDt = p.minDt;
      iam = p.iam;
      tag = p.tag;
      return *this;
   }
   
   bool is_neighbr(baseParticle &_p) {
      return (&_p != this) && (len2(r-_p.r) <= pow(KERNAL_RADIUS*max(h,_p.h),2));
   }

   bool isScaledNeighbr(baseParticle &_p,const double scale) {
      return (&_p != this) && (len2(r-_p.r) <= pow(KERNAL_RADIUS*scale*max(h,_p.h),2));
   }
   double mass;
   vect r;
   vect v;
   vect vhat;
   vect f;
   double h;
   double dens;
   double minDt;
   enum iamTypes iam;
   int tag;

   RETVARIABLE(position,vect,baseParticle,p.r)
   RETVARIABLE(velocity,vect,baseParticle,p.v)
   RETVARIABLE(vhat,vect,baseParticle,p.vhat)
   RETVARIABLE(density,double,baseParticle,p.dens)
   RETVARIABLE(smoothingLength,double,baseParticle,p.h)
   RETVARIABLE(tag,int,baseParticle,p.tag)
   RETVARIABLE(iam,iamTypes,baseParticle,p.iam)
   RETVARIABLE(mass,double,baseParticle,p.mass)
   RETVARIABLE(totalForce,vect,baseParticle,p.f)
};

inline baseGhost& baseGhost::operator=(const baseParticle &p) {
         r = p.r;
         v = p.v;
         vhat = p.vhat;
         mass = p.mass;
         dens = p.dens;
         h = p.h;
         tag = p.tag;
         iam = p.iam;
         return *this;
};

inline baseParticle& baseParticle::operator=(const baseGhost &g) {
      r = g.r;
      v = g.v;
      vhat = g.vhat;
      h = g.h;
      mass = g.mass;
      dens = g.dens;
      tag = g.tag;
      iam = g.iam;
      return *this;
};



#endif
