#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

inline void addGravity(Cparticle &p,CglobalVars &g) {
   p.f[1] -= 9.81;
}

void CcustomSim::beforeMiddle(double newTime) { 
   if (data->globals.sphStep > 400) {
      data->traverse<addGravity,Nsph::ifSphOrDem>();
   } else {
      data->traverse<addGravity,Nsph::ifSph>();
   }

}

