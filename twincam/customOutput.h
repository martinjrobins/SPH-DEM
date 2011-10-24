#ifndef CUSTOMOUTPUT_H
#define CUSTOMOUTPUT_H

#include <gsl/gsl_rng.h>
#include <ext/slist>
#include "customOutputBase.h"
#include "customSim.h"
#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"

namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later

class CpairData {
   public:
      Cparticle *p1,*p2;
      double initR2,initT,finalR2;
      ofstream *outFile;
};

class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {
         //_data->globals.custom.resize(2);

         rng = gsl_rng_alloc(gsl_rng_ranlxd1);
         gsl_rng_set(rng,time(NULL));

      }
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);
   private:
      //pair doubling time data
      Sgi::slist<CpairData> pairData;
      vector<ofstream *> openFiles;
      gsl_rng *rng;


};


#endif
