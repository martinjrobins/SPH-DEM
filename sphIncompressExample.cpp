#include "sphIncompressExample.h"

void CsphIncompress::start() {
   data->traverse<driftRAndKick,ifSphOrSphBoundaryOrDem>();
   data->traverse<driftRest,ifSphOrSphBoundary>();
   data->globals.time += data->globals.dt/2;

void CsphIncompress::middle() {
   data->globals.newDt = data->globals.maxdt;

   data->reset();
   data->traverse<calcPressSpsoundPdr2,ifSphOrSphBoundaryOrGhost>();	

   timeval t1,t2;
   gettimeofday(&t1,NULL);
   data->neighbours<calcForce,ifSphOrSphBoundary>();
   data->traverse<accelCondition,ifSphOrDem>(); 
   gettimeofday(&t2,NULL);
   data->globals.wtTotalForceCalc.tv_sec += t2.tv_sec-t1.tv_sec;
   data->globals.wtTotalForceCalc.tv_usec += t2.tv_usec-t1.tv_usec;

   double oldDt = data->globals.dt;
   if (data->globals.dt == 0.0) {
      if (data->globals.mpiRank==0) cout << "I think that this is the first timestep, reinitiallizing dt...."<<endl;
      data->globals.dt = data->globals.newDt;
   } else {
      data->globals.dt = data->globals.newDt;
   }
   data->setGlobalTimestep(data->globals.dt);
   data->traverse<kick,ifSphOrSphBoundaryOrDem>(0.5*oldDt,0.5*data->globals.dt);
}

void CsphIncompress::end() {
   data->traverse<driftR,ifSphOrSphBoundaryOrDem>();
   data->globals.time += data->globals.dt/2;
   data->reset();
   timeval t1,t2;
   gettimeofday(&t1,NULL);
   data->neighbours<calcDddtDudt,ifSphOrSphBoundary>();
   gettimeofday(&t2,NULL);
   data->globals.wtTotalDddtCalc.tv_sec += t2.tv_sec-t1.tv_sec;
   data->globals.wtTotalDddtCalc.tv_usec += t2.tv_usec-t1.tv_usec;
   data->traverse<driftRest,ifSphOrSphBoundary>();
   data->traverse<calcEnergies,ifSphOrSphBoundary>();
}
