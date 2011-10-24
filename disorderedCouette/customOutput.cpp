#include "customOutput.h"
#include "sph.h"
#include "misc.impl.h"
#include "dataLL.impl.h"

using namespace Nsph;

inline void calcEnstrophyAndAveDens(Cparticle &p,CglobalVars &g) {
   g.custom[0] += 0.5*p.mass*pow(p.vort,2)/p.dens;
   g.custom[1] += p.dens;
}


void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {

#ifdef VORT_LEASTSQUARES
   data->neighboursGroup<calcVortLeastSquares,ifSph>();
#else
   data->neighboursGroup<calcVortSPHSum,ifSph>();
#endif
   data->globals.custom[0] = 0;
   data->globals.custom[1] = 0;
   data->traverse<calcEnstrophyAndAveDens,ifSphOrSphBoundary>();
   data->globals.custom[1] /= data->numParticles();
}


inline bool ifInteriorSph(Cparticle &p) {
   return (p.iam == sph)&&(p.r[1]<RMIN[1]+0.8*(RMAX[1]-RMIN[1]))&&(p.r[1]>RMIN[1]+0.2*(RMAX[1]-RMIN[1]));
}

void CcustomOutput::calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io) {

   calcVelocityStructure<ifInteriorSph>(data,rng,io);
   
   calcPositionStructure<ifInteriorSph>(data,rng,io);

   outputVelocitySpectrum<ifSph>(data,io);
   
}
