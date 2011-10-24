#include "customOutput.h"
#include "sph.h"
#include "misc.impl.h"
#include "dataLL.impl.h"

using namespace Nsph;

void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {
#ifdef VORT_LEASTSQUARES
   data->neighboursGroup<calcVortLeastSquares,ifSph>();
#else
   data->neighboursGroup<calcVortSPHSum,ifSph>();
#endif
}

bool ifNotMiddleParticle(Cparticle &p) {
   return (p.tag!=-1);
}

void CcustomOutput::calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io) {
          
   outputVelocitySpectrum<ifNotMiddleParticle>(data,io);
   
}
