#ifndef CUSTOMSIM_H
#define CUSTOMSIM_H

#include "customSimBase.h"
#include "customConstants.h"
#include "dataLL.h"
#include "sph.h"
#include "vect.h"


inline void setDensity(Cparticle &p,CglobalVars &g) {
   p.dens = p.porosity*p.dens;
   p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
}
class CcustomSim : public CcustomSimBase {
   public:
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
         data->neighboursGroup<Nsph::calcPorosity,Nsph::ifSph>();
         data->traverse<setDensity,Nsph::ifSph>();
      }
      
      void beforeMiddle(double newTime);
};


#endif
