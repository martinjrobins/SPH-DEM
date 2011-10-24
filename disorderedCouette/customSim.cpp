#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

inline void driftR(Cparticle &p,CglobalVars &g,double dt) {
   p.r += dt*p.v;
}

void applyForcing(Cparticle &p,CglobalVars &g,vector<CfMode> &forceModes) {
   p.ff = 0;
   const double mult = 2.0*PI/(RMAX[1]-RMIN[1]);
   for (vector<CfMode>::iterator i=forceModes.begin();i!=forceModes.end();i++) {
      p.ff[0] += i->c*i->k[1]*sin(mult*i->k[0]*(p.r[0]-RMIN[0])+ mult*i->k[1]*(p.r[1]-RMIN[1]) + i->p);
      p.ff[1] += -i->c*i->k[0]*sin(mult*i->k[0]*(p.r[0]-RMIN[0])+ mult*i->k[1]*(p.r[1]-RMIN[1])+ i->p);
   }
   const double fmult = (1-exp(-100.0*pow(1.0-pow(p.r[0]/RMAX[0],2),2))) * (1-exp(-100.0*pow(1.0-pow(p.r[1]/RMAX[1],2),2)));
   p.ff *= fmult;

   p.f[0] += p.ff[0];
   p.f[1] += p.ff[1];
}


void CcustomSim::setupForcing() {
   vectInt k;
   for (k[0]=-FKMAX;k[0]<=FKMAX;k[0]++) {
      for (k[1]=-FKMAX;k[1]<=FKMAX;k[1]++) {
         double klen2 = len2(k);
         double klen = sqrt(klen2);
         if ((klen>=FKMIN)&&(klen<=FKMAX)&&(k[0]!=0)&&(k[1]!=0)) {
            CfMode newfMode;
            newfMode.k = k;
            newfMode.c = 0;
            newfMode.p = 0;
            forceModes.push_back(newfMode);
         }
      }
   }
   rng = gsl_rng_alloc(gsl_rng_ranlxd1);
   gsl_rng_set(rng,1);
   //gsl_rng_set(rng,time(NULL));
}

void CcustomSim::iterateForcing() {
   const double tmp = data->globals.dt/(2.0*FORCE_TAU);
   const double sigma = (1.0-tmp)/(1.0+tmp);
   for (vector<CfMode>::iterator i=forceModes.begin();i!=forceModes.end();i++) {
      const double rand = gsl_ran_gaussian(rng,1);
      const double amp = FORCE_AMP*(RMAX[1]-RMIN[1])/(2.0*PI*len2(i->k));
      double realPart = i->c*cos(i->p);
      double imagPart = i->c*sin(i->p);
      realPart = sqrt(1.0-pow(sigma,2))*amp*cos(PI*rand) + sigma*realPart;
      imagPart = sqrt(1.0-pow(sigma,2))*amp*sin(PI*rand) + sigma*imagPart;
      i->p = atan2(imagPart,realPart);
      i->c = sqrt(pow(realPart,2)+pow(imagPart,2));
      //if ((i->k[0]==8)&&(i->k[1]==1)) {
      //   i->c = 1;
      //} else {
      //   i->c = 0;
      //}
   }
}

//void CcustomSim::start(double newTime) {

//}

void CcustomSim::beforeMiddle(double newTime) {
   if (newTime<MAXTIME/OUTSTEP*FORCE_TILL) {
      iterateForcing();
      data->traverse<vector<CfMode>,applyForcing,Nsph::ifSph>(forceModes);
   }

   double dt = newTime-time;
   time = newTime;
   data->traverse<driftR,Nsph::ifSphBoundary>(dt);
}

void CcustomSim::beforeEnd(double newTime) {
   double dt = newTime-time;
   time = newTime;
   data->traverse<driftR,Nsph::ifSphBoundary>(dt);
}
