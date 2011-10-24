#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"
#include "particle.h"

inline void driftR(Cparticle &p,CglobalVars &g,double dt) {
   p.r += dt*p.v;
}

void CcustomSim::beforeMiddle(double newTime) {
   double dt = newTime-time;
   time = newTime;
   data->traverse<driftR,Nsph::ifSphBoundary>(dt);
}
void CcustomSim::beforeEnd(double newTime) {
   double dt = newTime-time;
   time = newTime;
   data->traverse<driftR,Nsph::ifSphBoundary>(dt);
}
