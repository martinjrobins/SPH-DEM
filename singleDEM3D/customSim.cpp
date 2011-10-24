#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

inline void addGravity(Cparticle &p,CglobalVars &g) {
   p.f[2] -= 9.81;
}

void CcustomSim::beforeMiddle(double newTime) { 
   //if (data->globals.sphStep > DAMP_STEPS) {
   //   data->traverse<addGravity,Nsph::ifDem>();
   //}
   data->traverse<addGravity,Nsph::ifSphOrDem>();
}

