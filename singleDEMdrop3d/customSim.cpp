#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

inline void addGravity(Cparticle &p,CglobalVars &g) {
   p.f[2] -= 9.81;
}

void CcustomSim::beforeMiddle(double newTime) { 
   if (data->globals.sphStep > DAMP_STEPS-100) {
      data->traverse<addGravity,Nsph::ifSphOrDem>();
   } else {
      data->traverse<addGravity,Nsph::ifSph>();
   }

}

