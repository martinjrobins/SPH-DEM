#include "customOutput.h"
#include "sph.h"
#include "dataLL.impl.h"
#include <fstream>

ofstream CcustomOutput::fo;

void outputParticle(Cparticle &p,CglobalVars &g,double newTime) {
   CcustomOutput::fo << p.r[1] <<' '<<p.v[0]<<'\n';
}

void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io) {
   double newTime = data->globals.time;
   if (newTime >= outputTimes[outputTimesIndex]) {
      cout << "\tOutputting Couette velocity profile for Time = "<<newTime;
      fo << "# Time = " << newTime << '\n';
      fo << "#   y   v_x\n";
      foExact << "# Time = " << newTime << '\n';
      foExact << "#   y   v_x\n";

      double ystart = RMIN[1]+2*PSEP;
      double ystop = RMAX[1]-2*PSEP;
      double yspan = ystop-ystart;

      for (double y=0;y<=yspan;y+=yspan/36.0) {
         double vx = VREF*y/yspan;
         for (int i=1;i<=20;i++) {
            vx += (2.0*VREF/(i*PI))*pow(-1.0,i)*sin(i*PI*y/(yspan))*exp(-VISCOSITY*pow(float(i),2)*pow(PI,2)*newTime/pow(yspan,2));
         }
         foExact << y+ystart << ' ' << vx << '\n';
      }
      foExact << '\n';
         
      data->traverse<outputParticle,Nsph::ifSphOrSphBoundary>(newTime);
      fo << '\n';
      outputTimesIndex++;
   }
}
