
#include <mpi.h>
#include "customConstants.h"
#include "sphIncompress.h"
#include "vect.h"
#include "dataLL.impl.h"
#include "io_data_vtk.h"
#include "particle.h"
#include <cstdlib>
#include <string>
#include <sys/time.h>
#include "customSim.h"
#include "customOutput.h"
#include "globalVars.h"

int main(int argc, char *argv[]) {

   if (argc != 4) {
      cout << "Usage: run infilename intimestep outfilename" << endl;
      return(-1);
   }
   cout << "Ran program with argv = ";
   for (int i=0;i<argc;i++) {
      cout << argv[i];
   }
   cout << endl;

   CglobalVars globals;
   particleContainer ps;

   //get rank and init all MPI stuff. Put rank in globals
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &(globals.mpiSize));
   MPI_Comm_rank(MPI_COMM_WORLD, &(globals.mpiRank));

   if (globals.mpiRank==0) cout << "Reading input file..." <<endl;
   string filename = argv[1];
   Cio_data_vtk io_data(filename,&globals);
   io_data.readDomain(atoi(argv[2]),&globals);
   io_data.readGlobals(atoi(argv[2]),&globals);
   io_data.readRestart(atoi(argv[2]),&ps,&globals);

   
   if (globals.mpiRank==0) cout << "after read, n = "<<ps.size()<<endl;
   
   int outstep = atoi(argv[2]);

   if (TIME_OVERIDE) {
      if (globals.mpiRank==0) cout << "time was set to time = "<<globals.time<<endl;
      globals.time = START_TIME;
      if (globals.mpiRank==0) cout << "Overiding time, time = "<<globals.time<<endl;
   }
   if (NSTEP_OVERIDE) {
      outstep = START_NSTEP;
      if (globals.mpiRank==0) cout << "Overiding output timestep  = "<<outstep<<endl;
   }
   
   string outFilename = argv[3];
   if (globals.mpiRank==0) cout << "Setting output filename as: "<<outFilename<<endl;
   io_data.setFilename(outFilename,&globals);


   if (globals.mpiRank==0) cout << "1time was set to time = "<<globals.time<<endl;

   if (globals.mpiRank==0) cout << "Creating data structure..." <<endl;

   /*
   CdataLL<sphParticle> *data = new CdataLL<sphParticle>(ps,globals,true);
#ifdef LIQ_DEM
   CdataLL<demParticle> *dataDem = new CdataLL<demParticle>(ps,globals,true);
   CsphLiqDem sph(data,dataDem);
   CcustomSim customSim(data,dataDem,globals.time);
   CcustomOutput customOutput(data,dataDem);
#else
   CsphIncompress sph(data);
   CcustomSim customSim(data,globals.time);
   CcustomOutput customOutput(data);
#endif
*/
   CdataLL *data = new CdataLL(ps,globals,true);
   CsphIncompress sph(data);
   CcustomSim customSim(data,globals.time);
   CcustomOutput customOutput(data);

   if (globals.mpiRank==0) cout << "2time was set to time = "<<globals.time<<endl;

   globals.maxdt = MAXTIME/OUTSTEP;

   if (globals.mpiRank==0) cout << "Starting simulation..."<<endl;
#ifndef START_WRITE_TIME
   double nextWriteTime = globals.time+MAXTIME/OUTSTEP;
#else
   double nextWriteTime = START_WRITE_TIME
#endif

   timeval wtStartSphStep,wtStartOutStep,wtStartSim,wtCurrTime;
   gettimeofday(&wtStartSim,NULL);
   wtStartOutStep = wtStartSim;
   wtStartSphStep = wtStartSim;
   globals.startOutStep();

   while ((globals.time <= MAXTIME) && (outstep <= OUTSTEP)) {

      globals.sphStep++;
      

      gettimeofday(&wtCurrTime,NULL);
      double simTimeSPHStep = wtCurrTime.tv_sec-wtStartSim.tv_sec+(wtCurrTime.tv_usec-wtStartSim.tv_usec)/1000000.0;
      //if (globals.mpiRank==0) cout << "\r#"<<globals.sphStep<<". Time = "<<globals.time<<". Simulation "<<globals.time/MAXTIME*100<<"\% complete. "<<wtCurrTime.tv_sec-wtStartSphStep.tv_sec+(wtCurrTime.tv_usec-wtStartSphStep.tv_usec)/1000000.0<<" sec per SPH step. Time taken = "<<simTimeSPHStep/(60*60)<<flush;
      wtStartSphStep = wtCurrTime;

      customSim.beforeStart(globals.time);
      sph.start();
      customSim.beforeMiddle(globals.time);
      sph.middle();
      customSim.beforeEnd(globals.time);
      sph.end();
      customSim.afterEnd(globals.time);
      
      if (globals.time>= nextWriteTime) {
      //if(globals.sphStep%10==0){
      //if (1) {
         cout << endl;
         outstep++;
         globals.outstep = outstep;
         globals.endOutStep();
         

         if (globals.mpiRank==0) cout << "\tWriting output at time = "<<globals.time<<". Simulation is "<<globals.time/MAXTIME*100<<"\% complete"<<endl;

         sph.calcOutputVars();
         customOutput.calcOutput(outstep,&customSim,&io_data);
         io_data.writeOutput(outstep,ps,&globals);
         io_data.writeGlobals(outstep,&globals);

         if (outstep%RESTART_EVERY_N_STEPS==0) {
            io_data.writeRestart(outstep,ps,&globals);
            io_data.writeDomain(outstep,&globals);
         }

         data->calcTraverseOrder();    //arrange traverse order spatially

         nextWriteTime = globals.time+MAXTIME/OUTSTEP;
         globals.maxdt = MAXTIME/OUTSTEP;


         gettimeofday(&wtCurrTime,NULL);
         double oSTime = wtCurrTime.tv_sec-wtStartOutStep.tv_sec+(wtCurrTime.tv_usec-wtStartOutStep.tv_usec)/1000000.0;
         double simTime = wtCurrTime.tv_sec-wtStartSim.tv_sec+(wtCurrTime.tv_usec-wtStartSim.tv_usec)/1000000.0;
         if (globals.mpiRank==0) cout <<'\t'<<oSTime<< " seconds per output timestep, "<<oSTime*(OUTSTEP-outstep)/(60.0*60.0)<<" hours remaining. "<<simTime/(60.0*60.0)<<" hours since start"<<endl;
         wtStartOutStep = wtCurrTime;

         if (globals.mpiRank==0) cout <<'\t'<<"Step profile: Total time = "<<globals.wtTotalOutStep.tv_sec+globals.wtTotalOutStep.tv_usec/1000000.0<<" Force = "<<globals.wtTotalForceCalc.tv_sec+globals.wtTotalForceCalc.tv_usec/1000000.0<<" dddt = "<<globals.wtTotalDddtCalc.tv_sec+globals.wtTotalDddtCalc.tv_usec/1000000.0<<" reset = "<<globals.wtTotalReset.tv_sec+globals.wtTotalReset.tv_usec/1000000.0<<" update domain = "<<globals.wtTotalUpdateDomain.tv_sec+globals.wtTotalUpdateDomain.tv_usec/1000000.0<<" MPI = "<<globals.wtTotalMPI.tv_sec+globals.wtTotalMPI.tv_usec/1000000.0;
#ifdef SMOOTHING
         if (globals.mpiRank==0) cout <<" Smooth = "<<globals.wtTotalSmooth.tv_sec+globals.wtTotalSmooth.tv_usec/1000000.0<<" SmoothMPI = "<<globals.wtTotalSmoothMPI.tv_sec+globals.wtTotalSmoothMPI.tv_usec/1000000.0<<endl;
#else
         if (globals.mpiRank==0) cout <<endl;
#endif
         globals.startOutStep();
         //exit(-1);
      }
   }
   MPI_Finalize();
}

