#ifndef CUSTOMSIM_H
#define CUSTOMSIM_H

#include "customSimBase.h"
#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define FKMIN 8
#define FKMAX 9

class CfMode {
   public:
      vectInt k;
      double c;
      double p;
};

class CcustomSim : public CcustomSimBase {
   public:
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
         setupForcing();
      }
      //void start(double newTime) {}
      void beforeMiddle(double newTime);
      void beforeEnd(double newTime);
      void setupForcing();
      void iterateForcing();
   private:
      vector<CfMode> forceModes;
      gsl_rng *rng;
};


#endif
