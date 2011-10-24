#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"


inline void setDensity(Cparticle &p,CglobalVars &g) {
   p.dens = p.r[1];
}


void CcustomSim::beforeMiddle(double newTime) {
}

void CcustomSim::beforeStart(double newTime) {
}

