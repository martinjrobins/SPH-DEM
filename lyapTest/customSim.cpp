#include "customSim.h"
#include "sph.h"
#include "dataLL.h"

inline void setVelocity(Cparticle &p,CglobalVars &g) {
   p.v = -VREF*cos(2.0*PI*p.r[0])*sin(2.0*PI*p.r[1]),VREF*sin(2.0*PI*p.r[0])*cos(2.0*PI*p.r[1]);
   p.v0 = p.v;
   p.vhat = p.v;
   p.vhat0 = p.v;
}


void CcustomSim::beforeMiddle(double newTime) { 
   data->traverse<setVelocity,Nsph::ifSph>();
}

