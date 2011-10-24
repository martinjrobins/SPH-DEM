#ifndef CUSTOMSIM_H
#define CUSTOMSIM_H

#include "customSimBase.h"
#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define FKMIN 7
#define FKMAX 9

class CfMode {
   public:
      vectInt k;
      double c;
};



class CcustomSim : public CcustomSimBase {
   public:
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
         frameAng = 0.0;
      }
      double getFrameAng() { return frameAng; }
      void beforeMiddle(double newTime);
      void beforeEnd(double newTime);
      int leftCamN;
      int rightCamN;
      int camN;
   private:
      void updateSim(const double newTime);
      double frameAng;
};


#endif
