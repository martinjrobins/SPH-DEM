#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

void findMiddleParticle(Cparticle &p,CglobalVars &g,int &pTag) {
   vect middle = 0.0;
   if ((pTag==-1)&&(len(p.r-middle)<=sqrt(0.5)*PSEP)) {
      pTag = p.tag;
      p.tag = -1;
      p.v = -JIG_A*2*PI*JIG_FREQ;
   }
}
void CcustomSim::initThisSim() {
   pTag = -1;
   data->traverse<int,findMiddleParticle,Nsph::ifSph>(pTag);
   if (pTag == -1) {
      cerr << "Error in CcustomSim: cannot find middle particle!"<<endl;
      exit(-1);
   }
}
bool ifMiddleParticle(Cparticle &p) {
   return (p.tag==-1);
}

void setMiddleParticle(Cparticle &p,CglobalVars &g,double time) {
   p.ff = -JIG_A*pow(2*PI*JIG_FREQ,2)*cos(2*PI*JIG_FREQ*time+PI/2);
   p.f += p.ff;
}

void CcustomSim::beforeMiddle(double newTime) {
   data->traverse<setMiddleParticle,ifMiddleParticle>(newTime);
   time = newTime;
}

