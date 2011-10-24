#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

inline void setVelocity(Cparticle &p,CglobalVars &g) {
   p.v = -VREF*cos(2.0*PI*p.r[0])*sin(2.0*PI*p.r[1]),VREF*sin(2.0*PI*p.r[0])*cos(2.0*PI*p.r[1]);
   p.v0 = p.v;
   p.vhat = p.v;
   p.vhat0 = p.v;
}


void CcustomSim::beforeMiddle(double newTime) { 
   if ((newTime<MAXTIME/OUTSTEP*20)&&(newTime+data->globals.dt >= MAXTIME/OUTSTEP*20)) {
      data->traverse<setVelocity,Nsph::ifSph>();
   }
   if ((newTime<MAXTIME/OUTSTEP*40)&&(newTime+data->globals.dt >= MAXTIME/OUTSTEP*40)) {
      data->traverse<setVelocity,Nsph::ifSph>();
   }
}

