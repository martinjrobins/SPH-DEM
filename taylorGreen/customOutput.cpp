#include "customOutput.h"
#include "dataLL.impl.h"
#include "sph.h"
#include "misc.h"

using namespace Nsph;


inline void calcMaxV(Cparticle &p,CglobalVars &g) {
   double v = len(p.v);
   double vhat = len(p.vhat);
   if (g.custom[0]<v) g.custom[0]=v;
   if (g.custom[1]<vhat) g.custom[1]=vhat;
}
inline void setColourToExactPress(Cparticle &p,CglobalVars &g) {
   p.colour = -0.25*(cos(4.0*PI*p.r[0])+cos(4.0*PI*p.r[1]))*exp(-4.0*4.0*PI*PI*g.time/REYNOLDS_NUMBER);
}


void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {
   data->neighboursGroup<calcVortLeastSquares,ifSph>();
   data->traverse<setColourToExactPress,ifSph>();
   
   data->globals.custom[0] = 0;
   data->globals.custom[1] = 0;
   data->traverse<calcMaxV,ifSphOrSphBoundary>();
}





     
