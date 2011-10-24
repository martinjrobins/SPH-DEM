
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
#include "misc.h"
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

   //get rank and init all MPI stuff. Put rank in globals
   //argc -= 3;
   //argv += 3;



   CglobalVars globals1,globals2;
   Cio_data_vtk io_data1(filename,&globals1);
   Nmisc::setupEntireDomain(globals1.procDomain,globals1.procNeighbrs);
   Cio_data_vtk io_data2(filename,&globals2);
   Nmisc::setupEntireDomain(globals2.procDomain,globals2.procNeighbrs);
   particleContainer ps1,ps2;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &(globals1.mpiSize));
   MPI_Comm_rank(MPI_COMM_WORLD, &(globals1.mpiRank));

   cout <<"reading 1st glboals"<<endl;
   io_data1.readGlobals(startStep,&globals1);
   globals1.outstep = startStep;
   cout <<"reading 1st outopu"<<endl;
   io_data1.readOutput(startStep,&ps1,&globals1);
   cout <<"creating data structure"<<endl;
   CdataLL *data1 = new CdataLL(ps1,globals1,true);
   data1->traverse<setMass>();

   io_data2.readGlobals(endStep,&globals2);
   globals2.outstep = endStep;
   io_data2.readOutput(endStep,&ps2,&globals2);
   CdataLL *data2 = new CdataLL(ps2,globals2,true);
   data2->traverse<setMass>();

   calcPoincare(data1,data2,&io_data1,20,1000);
}

