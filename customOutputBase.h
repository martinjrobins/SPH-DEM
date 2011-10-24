#ifndef CUSTOMOUTPUTBASE_H
#define CUSTOMOUTPUTBASE_H

#include "customConstants.h"
#include "customSim.h"
#include "io_data_vtk.h"
#include "dataLL.h"
#include "vect.h"

class CcustomOutputBase {
   public:
      CcustomOutputBase(CdataLL *_data): data(_data) {}
      virtual void calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {}
   protected:
      CdataLL *data;
};


#endif
