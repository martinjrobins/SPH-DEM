#ifndef CUSTOMOUTPUT_H
#define CUSTOMOUTPUT_H

#include "customOutputBase.h"
#include "customSim.h"
#include "customConstants.h"
#include "vect.h"


class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {
         foDem.open("demVel.dat");
         foSph.open("sphVel.dat");
         foDem << "time vx vy vz fp_x fp_y fp_z fv_x fv_y fv_z fb_x fb_y fb_z ff_x ff_y ff_z porosity rx ry rz fvelx fvely fvelz"<<endl;
         foSph << "time vx vy vz fp_x fp_y fp_z fv_x fv_y fv_z fb_x fb_y fb_z ff_x ff_y ff_z porosity rx ry rz fvelx fvely fvelz"<<endl;
      }
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);

      static ofstream foSph;
      static ofstream foDem;
   private:

};


#endif
