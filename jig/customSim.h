#ifndef CUSTOMSIM_H
#define CUSTOMSIM_H

#include "customSimBase.h"
#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


class CcustomSim : public CcustomSimBase {
   public:
      void initThisSim();
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
         initThisSim();
      }
      void beforeMiddle(double newTime);
   private:
      int pTag;
};


#endif
