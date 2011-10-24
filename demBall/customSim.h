#ifndef CUSTOMSIM_H
#define CUSTOMSIM_H

#include "customSimBase.h"
#include "customConstants.h"
#include "dataLL.h"
#include "sph.h"
#include "vect.h"


/*
inline void setDensity(Cparticle &p,CglobalVars &g) {
   p.dens = p.porosity*p.dens;
   p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
}
*/
class CcustomSim : public CcustomSimBase {
   public:
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
         //data->neighboursGroup<Nsph::calcPorosity,Nsph::ifSph>();
         //data->traverse<setDensity,Nsph::ifSph>();
#ifndef SETTLE
         data->globals.time = 0;
         Cparticle p;
         p.tag = data->getParticles()->size()+1;
         p.r[0] = (RMAX[0]-RMIN[0])/2;
         p.r[1] = (RMAX[1]-RMIN[1])/2;
         p.r[2] = RMAX[2];
         p.dens = DEM_DENS;
         p.h = DEM_RADIUS;
         p.mass = (4.0/3.0)*PI*pow(DEM_RADIUS,3)*DEM_DENS;
         p.v = 0.0;
         p.vhat = p.v;
         p.iam = dem;
         //data->insertNewParticle(p);
#endif
         first = true;
      }
      
      void beforeMiddle(double newTime);
      void beforeEnd(double newTime);
   private:
      bool first;
};


#endif
