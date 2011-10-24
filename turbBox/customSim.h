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
      double p;
};

class CcustomSim : public CcustomSimBase {
   public:
      CcustomSim(CdataLL *_data,double _time): CcustomSimBase(_data,_time) {
         frameAng = 0.0;
#ifdef FORCING
         setupForcing();
#endif
      }
      //void start(double newTime) {}
      double getFrameAng() { return frameAng; }
      void beforeMiddle(double newTime);
      void beforeEnd(double newTime);
      void setupForcing();
      void iterateForcing();
#ifdef REFERENCE_FRAME
      vect vFilter(const vect vin, const vect rin) { 
         vect vout;
         double angVel = calcAngVel(time);
         vout[0] = vin[0] + angVel*rin[1];
         vout[1] = vin[1] - angVel*rin[0];
         return vout;
      }
      vect rFilter(const vect rin) { 
         vect rout;
        double cang = cos(-frameAng);
         double sang = sin(-frameAng);
         rout[0] = rin[0]*cang-rin[1]*sang;
         rout[1] = rin[0]*sang+rin[1]*cang;
         return rout;
      }
#endif
   private:
      void updateSim(const double newTime);
      double calcAngVel(double time) {
         if (time<BOX_GLOB_TIME) {
            return BOX_GLOB_ANGVEL;
         } else {
            return BOX_GLOB_ANGVEL+BOX_PERT_ANGVEL*sin(2*PI*BOX_PERT_FREQ*(time-BOX_GLOB_TIME));
         }
      }
      double frameAng;
      vector<CfMode> forceModes;
      gsl_rng *rng;
      //double time;
      //CdataLL *data;
};


#endif
