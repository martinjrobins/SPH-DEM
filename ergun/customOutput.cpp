#include "customOutput.h"
#include "dataLL.impl.h"
#include "sph.h"
#include "misc.h"

using namespace Nsph;


ofstream CcustomOutput::foDem;
ofstream CcustomOutput::foSph;
inline void printVelSph(Cparticle &p,CglobalVars &g) {
   CcustomOutput::foSph << g.time <<' '<<p.vhat[0]<<' '<<p.vhat[1]<<' '<<p.vhat[2]<<' '<<p.fp[0]<<' '<<p.fp[1]<<' '<<p.fp[2]<<' '<<p.fv[0]<<' '<<p.fv[1]<<' '<<p.fv[2]<<' '<<p.fb[0]<<' '<<p.fb[1]<<' '<<p.fb[2]<<' '<<p.ff[0]<<' '<<p.ff[1]<<' '<<p.ff[2]<<' '<<p.porosity<<' '<<p.r[0]<<' '<<p.r[1]<<' '<<p.r[2]<<' '<<p.fvel[0]<<' '<<p.fvel[1]<<' '<<p.fvel[2]<<endl;
}
inline void printVelDem(Cparticle &p,CglobalVars &g) {
   CcustomOutput::foDem << g.time <<' '<<p.vhat[0]<<' '<<p.vhat[1]<<' '<<p.vhat[2]<<' '<<p.fp[0]<<' '<<p.fp[1]<<' '<<p.fp[2]<<' '<<p.fv[0]<<' '<<p.fv[1]<<' '<<p.fv[2]<<' '<<p.fb[0]<<' '<<p.fb[1]<<' '<<p.fb[2]<<' '<<p.ff[0]<<' '<<p.ff[1]<<' '<<p.ff[2]<<' '<<p.porosity<<' '<<p.r[0]<<' '<<p.r[1]<<' '<<p.r[2]<<' '<<p.fvel[0]<<' '<<p.fvel[1]<<' '<<p.fvel[2]<<endl;
}
void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {
   data->traverse<printVelDem,Nsph::ifDem>();
   data->traverse<printVelSph,Nsph::ifSph>();
}





     
