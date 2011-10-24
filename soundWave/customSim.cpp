#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"




#ifdef JIG

void damp(Cparticle &p,CglobalVars &g) {
   p.v *= 0.9;
}

const double sigma = 2*PSEP;
const double middle[2] = {(RMAX[0]-RMIN[0])/2.0+RMIN[0],(RMAX[1]-RMIN[1])/2.0+RMIN[1]};

void setVelocity(Cparticle &p,CglobalVars &g) {
   p.v[0] = 0.1*exp(-pow(p.r[0]-middle[0],2)/(2.0*pow(sigma,2)));
}

#endif

void CcustomSim::updateSim(const double newTime) {
   time = newTime;
}

void CcustomSim::densReit() {
   if (data->globals.mpiRank==0) cout<<"\tRe-initialising density..."<<endl;
   for (int i=0;i<10;i++) {
      data->reset();
      data->traverse<Nsph::calcHFromDens,Nsph::ifSphOrSphBoundary>();
      data->neighboursGroup<Nsph::calcDensity,Nsph::ifSphOrSphBoundary>();
   }
}

void CcustomSim::beforeMiddle(double newTime) {
   updateSim(newTime);
#ifdef JIG
   if (damping) {
      densReit();
   }
   if ((data->globals.time > DAMP_TIME)&&(damping)) {
      data->traverse<setVelocity,Nsph::ifSph>();   
      damping = false;
   }
#endif
}

void CcustomSim::beforeEnd(double newTime) {
#ifdef JIG
   if (damping) {
      data->traverse<damp,Nsph::ifSph>();
   }
#endif
   updateSim(newTime);
}

