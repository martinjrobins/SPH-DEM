#ifndef CUSTOMSIM_H
#define CUSTOMSIM_H

#include "customSimBase.h"
#include "customConstants.h"
#include "dataLL.h"
#include "sph.h"
#include "vect.h"


class CcustomSim : public CcustomSimBase {
   public:
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
         start = false;
         doOnce = false;
      }
      
      void beforeMiddle(double newTime);
      void beforeEnd(double newTime);
      void afterEnd(double newTime);
   private:
      bool start,doOnce;
};


#endif
