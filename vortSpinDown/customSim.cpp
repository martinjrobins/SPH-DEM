#include "customSim.h"
#include "sph.h"
#include "dataLL.h"

inline void setVelocity(Cparticle &p,CglobalVars &g) {
   p.v = 0.25*(p.r[1]-0.5),0.25*(0.5-p.r[0]);
}


void CcustomSim::beforeStart(double newTime) { 
   if ((newTime<MAXTIME/OUTSTEP*30)&&(newTime+data->globals.dt >= MAXTIME/OUTSTEP*30)) {
      data->traverse<setVelocity,Nsph::ifSph>();
   }
}

