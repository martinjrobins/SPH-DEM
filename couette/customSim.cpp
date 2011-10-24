#include "customSim.h"
#include "sph.h"
#include "particle.h"
#include "dataLL.impl.h"

inline void driftR(Cparticle &p,CglobalVars &g,double dt) {
   p.r += dt*p.v;
}

void CcustomSim::start(double newTime) {
   //time = newTime;
}
void CcustomSim::middle(double newTime) {
   double dt = newTime-time;
   time = newTime;
   data->traverse<driftR,Nsph::ifSphBoundary>(dt);
}
void CcustomSim::end(double newTime) {
   double dt = newTime-time;
   time = newTime;
   data->traverse<driftR,Nsph::ifSphBoundary>(dt);
}
