#ifndef CUSTOMOUTPUT_H
#define CUSTOMOUTPUT_H

#include "customOutputBase.h"
#include "customSim.h"
#include "customConstants.h"
#include "vect.h"


class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {
         fo.open("demVel.dat");
         fo << "time vx vy vz fp_x fp_y fp_z fv_x fv_y fv_z fb_x fb_y fb_z ff_x ff_y ff_z porosity rx ry rz"<<endl;
      }
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);

      static ofstream fo;
   private:

};


#endif
