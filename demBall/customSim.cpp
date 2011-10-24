#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

inline void addGravity(Cparticle &p,CglobalVars &g) {
   p.f[2] -= 9.81;
}
inline void addGravityAndParticleMass(Cparticle &p,CglobalVars &g) {
   p.f[2] -= 9.81;
   p.f[2] -= ((0.5*(RMAX[2]-RMIN[2])-PSEP)*(RMAX[0]-RMIN[0])*(RMAX[1]-RMIN[1]))*(1.0-POROSITY)*(DEM_DENS/DENS)*9.81;
}
inline void freezeParticles(Cparticle &p,CglobalVars &g) {
   p.f = 0.0;
   p.v = 0.0;
//#ifdef MANY_PARTICLES
//   p.v[2] = -VREF;
//#endif
   p.vhat = p.v;
}

inline void setFluidVel(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g) {
      double sum = 0;
      vect pvel = 0.0;
      int n = neighbrs.size();
      for (int i=0;i<n;i++) {
         Cparticle *pn = neighbrs[i];
         if (pn->iam!=dem) continue;
         const vect dx = p.r-pn->r;
         const double r = len(dx);
         const double dvol = pn->mass/pn->dens;
         const double q = r/pn->h;
         const double Wab = Nsph::W(q,pn->h);
         const double dvWab = dvol*Wab;
         pvel += dvWab*pn->vhat;
         sum += dvWab;
      }
      if (sum>0.0) {
         double rsum = 1/sum;
         pvel *= rsum;
      }
      p.vhat = -(1-p.porosity)*pvel;
      p.v = p.vhat/p.porosity;
}


void CcustomSim::beforeMiddle(double newTime) { 
   if (data->globals.time > TIME_DROP_PARTICLE) {
      data->traverse<addGravity,Nsph::ifSphOrDem>();
   } else {
      data->traverse<addGravity,Nsph::ifSph>();
      data->traverse<freezeParticles,Nsph::ifDem>();
   }
}

void CcustomSim::beforeEnd(double newTime) { 
   if (data->globals.time <= TIME_DROP_PARTICLE) {
      data->traverse<freezeParticles,Nsph::ifDem>();
   }
}
