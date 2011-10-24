#ifndef CUSTOMSIM_H
#define CUSTOMSIM_H

#include "customSimBase.h"
#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"

class CcustomSim : public CcustomSimBase {
   public:
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
      }
};


#endif
