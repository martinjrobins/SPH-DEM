
#include "customConstants.h"
#include "customOutput.h"
#include "vect.h"
#include "sphIncompress.h"
#include "dataLL.impl.h"
#include "io_data_vtk.h"
#include "particle.h"
#include "sph.h"
#include "customSim.h"
#include <cstdlib>
#include <string>

int main(int argc, char *argv[]) {
   if (argc != 4) {
      cout << "Usage: post infilename startStep endStep" << endl;
      return(-1);
   }
   string filename = argv[1];
   int startStep = atoi(argv[2]);
   int endStep = atoi(argv[3]);
   CglobalVars globals;
   Cio_data_vtk io_data(filename,&globals);
   particleContainer ps;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &(globals.mpiSize));
   MPI_Comm_rank(MPI_COMM_WORLD, &(globals.mpiRank));

   Nmisc::setupEntireDomain(globals.procDomain,globals.procNeighbrs);

   cout << "Processor "<<globals.mpiRank<<" of "<<globals.mpiSize<<" has domain ";
   for (int i=0;i<NDIM*2;i++) {
      cout<<globals.procDomain[i]<<' ';
   }
   cout << " and neighbrs "<<globals.procNeighbrs<<endl;
  
   cout <<"read first data set"<<endl;
   io_data.readCSIRO(startStep,&ps,&globals);
   cout <<"creating data structure"<<endl;
   CdataLL *data = new CdataLL(ps,globals,false);
   cout <<"creating custome sutff"<<endl;
   CcustomSim *customSim = new CcustomSim(data,globals.time);
   

   for (int nstep=startStep;nstep<=endStep;nstep++) {
      globals.outstep = nstep;
      globals.time = (nstep-1)*(MAXTIME/OUTSTEP);
      //cout <<"reading glboals"<<endl;
      //io_data.readGlobals(nstep,&globals);
      cout <<"reading outopu"<<endl;
      io_data.readCSIRO(nstep,&ps,&globals);
      io_data.writeOutput(nstep,ps,*customSim,&globals);
      io_data.writeGlobals(nstep,&globals);
      
      /*
      customSim.beforeStart(globals.time);
      customSim.beforeMiddle(globals.time);
      customSim.beforeEnd(globals.time);
      customSim.afterEnd(globals.time);
      customOutput.calcOutput(outstep,&customSim,&io_data);
      */ 

      cout << "Timestep "<<nstep<<" Time "<<globals.time<<" Conversion "<<double(nstep-startStep)/double(endStep-startStep)*100<<"\% complete"<<endl;
      
   }
}

