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
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
         damping = true;
      }
      void beforeMiddle(double newTime);
      void beforeEnd(double newTime);
   private:
      void updateSim(const double newTime);
      void densReit();
      bool damping;
};


#endif
