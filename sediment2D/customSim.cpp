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
   p.vhat = 0.0;
}

/*
inline void addCustomBoundaries(Cparticle &p,CglobalVars &g) {
   const double overlap = 2.0*DEM_RADIUS-r;
   if (p.r[2] < RMIN[2]+DEM_RADIUS) {
      const double overlap = 2.0*DEM_RADIUS-r;

   const vect dx = p.r-pn->r;
   const double r = len(dx);
   const vect dv = p.vhat - pn->vhat;
   const vect normal = dx/r;
   if (overlap>0) {
      const double overlap_dot = -dot(dv,normal);
      const vect fNorm = demNormal(overlap,overlap_dot,normal)/p.mass;
      //const double reducedRadius = 0.5*DEM_RADIUS; 
      //p.dVortDt += dem_tangential(tdv[0],normal,reducedRadius,DEM_K_SLIDING,DEM_GAMMA_SLIDING,genFtest(DEM_MU_STATIC,DEM_K_STATIC,fNorm,overlap),genFtest(DEM_MU_SLIDING,DEM_K_SLIDING,fNorm,overlap),p.eta[0],
      //const vect fslid = dem_tangential(dv - normal*(dot(normal,dv)),   );
      p.fp += fNorm;
      p.f += fNorm;
   }
}
*/



void CcustomSim::beforeMiddle(double newTime) { 
   //if (data->globals.sphStep > DAMP_STEPS) {
   //   data->traverse<addGravity,Nsph::ifDem>();
   //}
   if (data->globals.time > TIME_DROP_PARTICLE) {
      data->traverse<addGravity,Nsph::ifSphOrDem>();
   } else {
//#ifdef MANY_PARTICLES
   //   data->traverse<addGravityAndParticleMass,Nsph::ifSph>();
//#else
      data->traverse<addGravity,Nsph::ifSph>();
//#endif
   }
   //data->traverse<addCustomBoundaries,Nsph::ifDem>();
}

void CcustomSim::beforeEnd(double newTime) { 
   if (data->globals.time <= TIME_DROP_PARTICLE) {
      data->traverse<freezeParticles,Nsph::ifDem>();
   }
}
