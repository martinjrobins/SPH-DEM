#ifndef CUSTOMOUTPUT_H
#define CUSTOMOUTPUT_H

#include "customOutputBase.h"
#include "customSim.h"
#include "customConstants.h"
#include "vect.h"

class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {
         fo.open("com.dat");
         fo << "time dry0.3:r,h dry0.4 dry0.5 dry0.6 dry0.7 still10 still100 still1000 still10000 still100000"<<endl;
      }
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);

      ofstream fo;
};


#endif
