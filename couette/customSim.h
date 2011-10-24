#ifndef CUSTOMSIM_H
#define CUSTOMSIM_H


#include "customSimBase.h"
#include "customConstants.h"

class CcustomSim : public CcustomSimBase {
   public:
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {}
      void start(double newTime);
      void middle(double newTime);
      void end(double newTime);
};


#endif
