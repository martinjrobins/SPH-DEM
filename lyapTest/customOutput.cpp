#include "customOutput.h"
#include "sph.h"

using namespace Nsph;

inline void setMass(Cparticle &p,CglobalVars &g) {
   p.mass = PSEP*PSEP*DENS;
}

inline void calcMaxV(Cparticle &p,CglobalVars &g) {
   double v = len(p.v);
   if (g.custom[0]<v) g.custom[0]=v;
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

void CcustomOutput::calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io) {
   int oldstep = io->findClosestStepToTime(data->globals.time-LYAP_DT);
   //cout <<" current time is "<<data->globals.time<<endl;
   //cout <<" closest step to t = "<<data->globals.time-LYAP_DT<<" is "<<oldstep<<endl;
   if (oldstep >= 1) {
      CglobalVars g;
      g.outstep = oldstep;
      particleContainer ps;
      io->readGlobals(oldstep,&g);
      io->readOutput(oldstep,&ps,&g);
      Nmisc::setupEntireDomain(g.procDomain,g.procNeighbrs);

      cout <<"calculating lyap with dt = "<<data->globals.time-g.time<<endl;

      CdataLL oldData(ps,g);
      oldData.traverse<setMass>();

      Nmisc::calcLyap(&oldData,data,io);
   }
}
