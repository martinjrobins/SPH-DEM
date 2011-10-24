#include "customOutput.h"
#include "sph.h"
#include "dataLL.impl.h"

using namespace Nsph;

void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {
   data->neighboursGroup<calcVortLeastSquares,ifSph>();
}
