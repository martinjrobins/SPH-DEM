#ifndef CUSTOMOUTPUT_H
#define CUSTOMOUTPUT_H

#include "customOutputBase.h"
#include "customSim.h"
#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"

class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {}
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);
};


#endif
