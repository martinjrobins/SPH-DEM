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
         first = true;
      }
      
      void beforeStart(double newTime);
      void beforeMiddle(double newTime);
      void beforeEnd(double newTime);
      void afterEnd(double newTime);
   private:
      bool first;
};


#endif
