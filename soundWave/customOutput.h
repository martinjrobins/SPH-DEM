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
#include "sphIncompress.h"


using namespace Nmisc;

class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {
      }
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);
      void calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io, CsphIncompress *sph);
   private:
};


#endif
