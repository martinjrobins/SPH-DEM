#include "customSim.h"
#include "dataLL.impl.h"

const double RANDOMTAG = 955234;

inline void addGravity(Cparticle &p,CglobalVars &g) {
   p.ff[2] -= 9.81;
   p.f[2] -= 9.81;
}


inline void goWallgo(Cparticle &p,CglobalVars &g,double &dt) {
   if (p.r[0] > L1+L2+PSEP) return;
   p.v[2] = WALLSPEED;
   p.vhat[2] = WALLSPEED;
   p.r[2] += dt*WALLSPEED;
   if (p.r[2]>RMAX[2]) {
      p.tag = -111;
   }
}

void addParticle(Cparticle &p,CdataLL &data,double y,int tag) {
   Cparticle pnew=p;
   pnew.r[1] = y;
   pnew.tag = tag;
   data.insertNewParticle(pnew);
}

inline void myDriftR(Cparticle &p,CglobalVars &g) {
   const double dt = g.dt/2;
   p.r += dt*p.vhat;
}
inline void myReverseDriftR(Cparticle &p,CglobalVars &g) {
   const double dt = g.dt/2;
   p.r -= dt*p.vhat;
}

inline void addPeriodicBoundaries(Cparticle &p,CglobalVars &g,CdataLL &data) {
   if ((p.tag==RANDOMTAG)||(p.tag==-111)) return;
   if (g.mpiRank<=ENDPERIODICCPU) {
      const double maxY = ((NY/2.0)+PY+1)*PSEP+RMIN[1];
      const double minY = ((NY/2.0)-PY)*PSEP+RMIN[1];
      const double width = maxY-minY;
      const double expandX = NXP*PSEP;
      //cout <<"p.r = "<<p.r<<endl;
      if (p.tag==ENDTAG) {
         if (p.r[0]>expandX) return;
         if ((p.r[1]>maxY)||(p.r[1]<minY)) {
            data.markForDeletion(p);
            return;
         } else {
            p.tag = STARTTAG;
         }
      }
      if (p.r[1]>maxY) {
         p.r[1] = p.r[1]-width;
      } else if (p.r[1]<minY) {
         p.r[1] = p.r[1]+width;
      }
      if (p.r[0] > expandX-2*H) {
         for (int i=0;i<RY;i++) {
            if (p.r[0] > expandX) {
               p.tag = ENDTAG;
               addParticle(p,data,p.r[1]-(i+1)*width,ENDTAG);
               addParticle(p,data,p.r[1]+(i+1)*width,ENDTAG);
            } else {
               addParticle(p,data,p.r[1]-(i+1)*width,RANDOMTAG);
               addParticle(p,data,p.r[1]+(i+1)*width,RANDOMTAG);
            }
         }
      } else {
         if (p.r[1]>maxY-2.1*H) {
            addParticle(p,data,p.r[1]-width,RANDOMTAG);
         } 
         if (p.r[1]<minY+2.1*H) {
            addParticle(p,data,p.r[1]+width,RANDOMTAG);
         }
      }
   }
}

inline void removePeriodicBoundaries(Cparticle &p,CglobalVars &g,CdataLL &data) {
   if (p.tag==RANDOMTAG) {
      data.markForDeletion(p);
   }
}


/*
inline void createOneWayBoundary(Cparticle &p,CglobalVars &g) {
   if (g.mpiRank==ENDPERIODICCPU) {
      if (p.r[0] > NX*PSEP) {

   } else if (g.mpiRank==ENDPERIODICCPU+1) {
   }
}
*/

void CcustomSim::beforeMiddle(double newTime) {
   data->traverse<addGravity,Nsph::ifSph>();
   if ((RY>0)&&(data->globals.mpiRank<=ENDPERIODICCPU)) {
      data->traverse<CdataLL,addPeriodicBoundaries,Nsph::ifSphOrSphBoundaryOrBoundary>(*data);
      //data->traverse<CdataLL,createOneWayBoundary,Nsph::ifGhost>();
   }
   if (newTime > WALLUP) {
      double dt = newTime-oldTime;
      data->traverse<double,goWallgo,Nsph::ifBoundaryOrSphBoundary>(dt);
   }
   oldTime = newTime;
}
void CcustomSim::beforeEnd(double newTime) {
   data->traverse<addGravity,Nsph::ifSph>();
   if ((RY>0)&&(data->globals.mpiRank<=ENDPERIODICCPU)) {
      data->traverse<CdataLL,removePeriodicBoundaries,Nsph::ifSphOrSphBoundaryOrBoundary>(*data);
      data->traverse<myDriftR,Nsph::ifSphOrSphBoundary>();
      data->traverse<CdataLL,addPeriodicBoundaries,Nsph::ifSphOrSphBoundaryOrBoundary>(*data);
      data->traverse<myReverseDriftR,Nsph::ifSphOrSphBoundary>();
      //data->traverse<CdataLL,createOneWayBoundary,Nsph::ifGhost>();
   }
   if (newTime > WALLUP) {
      double dt = newTime-oldTime;
      data->traverse<double,goWallgo,Nsph::ifBoundaryOrSphBoundary>(dt);
   }
   oldTime = newTime;
}

void CcustomSim::afterEnd(double newTime) {
   data->traverse<addGravity,Nsph::ifSph>();
   if ((RY>0)&&(data->globals.mpiRank<=ENDPERIODICCPU)) {
      data->traverse<CdataLL,removePeriodicBoundaries,Nsph::ifSphOrSphBoundaryOrBoundary>(*data);
      data->deleteParticles();
   }
}
