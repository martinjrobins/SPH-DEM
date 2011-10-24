#ifndef CUSTOMOUTPUT_H
#define CUSTOMOUTPUT_H

//#include <fftw3.h>
#include <ext/slist>
#include "customOutputBase.h"
#include "customSim.h"
#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"
#include "misc.h"

namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later

using namespace Nmisc;

class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {

         rng = gsl_rng_alloc(gsl_rng_ranlxd1);
         gsl_rng_set(rng,time(NULL));

      }
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);
      void calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io);
   private:
      gsl_rng *rng;
};


#endif
