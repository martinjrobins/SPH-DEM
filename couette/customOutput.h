#ifndef CUSTOMOUTPUT_H
#define CUSTOMOUTPUT_H

#include "customOutputBase.h"
#include "customConstants.h"
#include <fstream>

class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {
         outputTimesIndex = 0;
         outputTimes[0] = 0.0100;
         outputTimes[1] = 0.0225;
         outputTimes[2] = 0.045;
         outputTimes[3] = 0.1125;
         outputTimes[4] = 0.350;
         outputTimes[5] = 2000;
         fo.open("couetteVelocityProfiles.dat");
         foExact.open("couetteVelocityProfilesExact.dat");
         fo << "# Times are about "<<outputTimes[0]<<' '<<outputTimes[1]<<' '<<outputTimes[2]<<' '<<outputTimes[3]<<' '<<outputTimes[4]<<". Real times will vary, see below.\n";
         foExact << "# Times are about "<<outputTimes[0]<<' '<<outputTimes[1]<<' '<<outputTimes[2]<<' '<<outputTimes[3]<<' '<<outputTimes[4]<<". Real times will vary, see below.\n";
      };
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);

      static ofstream fo;
   private:
      ofstream foExact;
      double outputTimes[6];
      int outputTimesIndex;
};


#endif
