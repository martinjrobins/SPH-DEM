
#include <mpi.h>
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

inline void setMass(Cparticle &p,CglobalVars &g) {
   p.mass = PSEP*PSEP*DENS;
}

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

   //get rank and init all MPI stuff. Put rank in globals
   //argc -= 3;
   //argv += 3;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &(globals.mpiSize));
   MPI_Comm_rank(MPI_COMM_WORLD, &(globals.mpiRank));

   Nmisc::setupEntireDomain(globals.procDomain,globals.procNeighbrs);

   cout << "Processor "<<globals.mpiRank<<" of "<<globals.mpiSize<<" has domain ";
   for (int i=0;i<NDIM*2;i++) {
      cout<<globals.procDomain[i]<<' ';
   }
   cout << " and neighbrs "<<globals.procNeighbrs<<endl;
  

   cout <<"reading 1st glboals"<<endl;
   io_data.readGlobals(startStep,&globals);
   cout <<"reading 1st outopu"<<endl;
   io_data.readOutput(startStep,&ps,&globals);
   cout <<"creating data structure"<<endl;
   CdataLL *data = new CdataLL(ps,globals,true);
   data->traverse<setMass>();
   cout <<"creating custome sutff"<<endl;
   CcustomSim *customSim = new CcustomSim(data,globals.time);
   CcustomOutput *customOutput = new CcustomOutput(data);
   CsphIncompress sph(data);

   for (int nstep=startStep;nstep<=endStep;nstep++) {
      globals.outstep = nstep;
      cout <<"reading glboals"<<endl;
      io_data.readGlobals(nstep,&globals);
      cout <<"reading outopu"<<endl;
      io_data.readOutput(nstep,&ps,&globals);
      data->traverse<setMass>();
      cout <<"updating particle"<<endl;
      data->updateParticles();
      data->reset();
      sph.calcOutputVars();
      
      /*
      customSim.beforeStart(globals.time);
      customSim.beforeMiddle(globals.time);
      customSim.beforeEnd(globals.time);
      customSim.afterEnd(globals.time);
      */ 
      customOutput->calcOutput(nstep,customSim,&io_data);

      cout << "Timestep "<<nstep<<" Time "<<globals.time<<" Post-Process "<<double(nstep-startStep)/double(endStep-startStep)*100<<"\% complete"<<endl;
      
      customOutput->calcPostProcess(nstep,customSim,&io_data,&sph);
      
   }
}

