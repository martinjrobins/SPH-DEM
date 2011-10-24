#include "customSim.h"
#include "dataLL.impl.h"

inline void addGravity(Cparticle &p,CglobalVars &g) {
   p.ff[2] -= 9.81;
   p.f[2] -= 9.81;
}

void CcustomSim::beforeMiddle(double newTime) {
   data->traverse<addGravity,Nsph::ifSph>();
}

