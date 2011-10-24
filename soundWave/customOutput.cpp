#include "customOutput.h"
#include "sph.h"
#include "misc.impl.h"
#include "dataLL.impl.h"

using namespace Nsph;


void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {

}

            
void CcustomOutput::calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io, CsphIncompress *sph) {
   vector<Cparticle> psGrid;
   vectInt gridDims;
   sph->renderToGrid(psGrid,gridDims);
   
   ofstream fo;
   string outputFN = "VGrid";
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   outputFN = io->getFilename()+outputFN+strTimestep+".dat";

   cout <<"writing grid data to file: "<<outputFN<<endl;
   fo.open(outputFN.c_str(),ios::out|ios::trunc);
   fo << "# "<<NX<<" x "<<NY<<" matrix"<<endl;

   for (particleContainer::iterator p = psGrid.begin();p != psGrid.end();p++) {
      if (p->r[0] == RMIN[0]) {
         fo<<endl;
      }
      fo<<p->v[0]<<' ';
   }
   fo.close();


   
   //outputVelocitySpectrum<ifSph>(data,io);
}
