#include "customOutput.h"
#include "sph.h"
#include "misc.impl.h"
#include "dataLL.impl.h"

using namespace Nsph;



void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {
}

            
void CcustomOutput::calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io, CsphIncompress *sph) {

          
   //calcVelocityStructure<ifInteriorSph>(data,rng,io);
   
   //calcPositionStructure<ifInteriorSph>(data,rng,io);

   //pairDoubling(data,pairData,openFiles,rng,io);
   
   outputVelocitySpectrum<ifSph>(data,io);
   
}
