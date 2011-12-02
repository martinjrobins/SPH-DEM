#include "sphIncompress.h"

void CsphIncompress::start() {
   data->traverse<driftRAndKick,ifSphOrSphBoundary>();
#ifdef LIQ_DEM
   data->traverse<driftR,ifDem>();
#endif


#ifdef SLK
   data->globals.slkIndex0++;
   data->globals.slkIndexLast++;
   if (data->globals.slkIndex0>=SLK_NUM_TIMESTEPS) data->globals.slkIndex0 = 0;
   if (data->globals.slkIndexLast>=SLK_NUM_TIMESTEPS) data->globals.slkIndexLast = 0;
#endif
   data->traverse<driftRest,ifSphOrSphBoundary>();
   data->globals.time += data->globals.dt/2;

   //if ((data->globals.sphStep < DAMP_STEPS)||(data->globals.sphStep%REINIT_DENS_EVERY_N_STEPS==0)) {
   if ((data->globals.sphStep%REINIT_DENS_EVERY_N_STEPS==0)) {
      data->reset();
#ifdef MY_VAR_RES
      data->neighboursGroup<calcGradVLeastSquares,ifSphOrSphBoundary>();
      data->syncParticlesBetweenProcs<vect,readGradV0,writeGradV0>();
      data->syncParticlesBetweenProcs<vect,readGradV1,writeGradV1>();
      //data->syncParticlesBetweenProcs<vect,readGradV2,writeGradV2>();
#endif
      if (data->globals.mpiRank==0) cout<<"\tRe-initialising density...sphStep="<<data->globals.sphStep<<endl;
      for (int i=0;i<10;i++) {
         data->traverse<calcHFromDens,ifSphOrSphBoundaryOrGhost>();
         data->syncParticlesBetweenProcs<double,readH,writeH>();
         data->neighboursGroup<calcDensity,ifSphOrSphBoundary>();
         data->syncParticlesBetweenProcs<double,readDens,writeDens>();
      }
   }

   data->traverse<initSumsMiddle,ifSphOrSphBoundary>();
}

void CsphIncompress::middle() {
   data->globals.newDt = data->globals.maxdt;

#ifdef SMOOTHED_VISC_VELOCITY
   data->reset(SMOOTH_HFAC_MULT);
   data->findBothNeighbours<ifSph>();
#endif

#ifdef SPH_SMOOTHING
   data->reset(SMOOTH_HFAC_MULT);
   data->findBothNeighbours<ifSph>();
#endif
#ifdef GRID_SMOOTHING
   data->reset();
#endif
#ifndef SMOOTHING
   //cout <<"rank = "<<data->globals.mpiRank<<"before reset "<<endl;
   data->reset();
   //cout <<"rank = "<<data->globals.mpiRank<<"after reset "<<endl;
#endif
#ifdef CORRECTED_GRADIENT
   data->neighboursGroup<calcInvM,ifSphOrSphBoundary>();
#endif
#ifdef VAR_H_CORRECTION 
   data->neighboursGroup<calcOmega,ifSphOrSphBoundary>();
   data->syncParticlesBetweenProcs<double,readOmega,writeOmega>();
#endif
   
   timeval t1,t2;
   gettimeofday(&t1,NULL);
#ifdef LIQ_DEM
   data->traverse<initPorosityAndDrag,ifSphOrSphBoundaryOrGhost>();
   data->neighboursGroupCoupling<calcPorosityAndDrag,ifDem>();
   data->reverseSyncParticlesBetweenProcs<vect,readFdrag,sumFdrag>();
   data->reverseSyncParticlesBetweenProcs<double,readPorosity,sumPorosity>();
   data->reverseSyncParticlesBetweenProcs<double,read_dPorositydt,sum_dPorositydt>();
   data->neighboursGroupAtRadius<finalisePorosityAndDrag,ifSphOrSphBoundary>(LIQ_DEM_COUPLING_RADIUS);
   data->traverse<setPorosityForBoundary,ifBoundary>();
   data->syncParticlesBetweenProcs<double,readPorosity,writePorosity>();
   /*
   data->neighboursGroupCoupling<calcPorosityAndDrag2,ifSphOrSphBoundary>();
   data->syncParticlesBetweenProcs<double,readPorosity,writePorosity>();
   */
#endif
   data->traverse<calcPressSpsoundPdr2,ifSphOrSphBoundaryOrGhost>();	

#ifdef SLK
   data->neighboursGroup<calcGLeastSquares,ifSphOrSphBoundary>();
   data->syncParticlesBetweenProcs<vectTensor,readG,writeG>();
#endif
#ifdef SMOOTHED_VISC_VELOCITY
   data->neighboursGroupScaled<smoothViscV,ifSphOrSphBoundary>();
   data->syncParticlesBetweenProcs<vect,readViscV,writeViscV>();
#endif
   //cout <<"rank = "<<data->globals.mpiRank<<"before force "<<endl;
   data->neighbours<calcForce,ifSphOrSphBoundary>();
   //cout <<"rank = "<<data->globals.mpiRank<<"after force"<<endl;
#ifdef BACKGROUND_PRESSURE_FORCE
   data->traverse<addBackgroundPressure,ifSphOrSphBoundary>();
#endif
#ifdef LIQ_DEM
   data->globals.newDt = min(data->globals.newDt,liqDemCondition());
#ifndef FIXED_DEM
#ifndef NO_DEM_CONTACTS
   if (data->globals.time > TIME_DROP_PARTICLE) {
      data->globals.newDt = min(data->globals.newDt,demCondition());
   }
#endif
#endif
#endif
   data->traverse<accelCondition,ifSphOrDem>(); 
   gettimeofday(&t2,NULL);
   data->globals.wtTotalForceCalc.tv_sec += t2.tv_sec-t1.tv_sec;
   data->globals.wtTotalForceCalc.tv_usec += t2.tv_usec-t1.tv_usec;

   double oldDt = data->globals.dt;
   if (data->globals.dt == 0.0) {
      if (data->globals.mpiRank==0) cout << "I think that this is the first timestep, reinitiallizing dt...."<<endl;
      data->globals.dt = data->globals.newDt;
   } else {
      //data->globals.dt = 2.0*data->globals.newDt - data->globals.dt;
      //data->globals.dt = data->globals.newDt*data->globals.dt/(2.0*data->globals.dt-data->globals.newDt);
      data->globals.dt = data->globals.newDt;
   }

   data->setGlobalTimestep(data->globals.dt);

   //data->neighboursGroup<calcVort,ifSphOrSphBoundary>();
   data->traverse<kick,ifSphOrSphBoundary>(0.5*oldDt,0.5*data->globals.dt);
   

   data->traverse<initSumsEnd,ifSphOrSphBoundary>();
#ifdef LIQ_DEM
   data->traverse<driftRAndKick,ifDem>();
   data->traverse<initSumsMiddle,ifDem>();
#endif
}

void CsphIncompress::end() {


#ifndef SMOOTHING
   data->traverse<driftR,ifSphOrSphBoundary>();
#endif
#ifdef SLK
   data->globals.slkIndex0++;
   data->globals.slkIndexLast++;
   if (data->globals.slkIndex0>=SLK_NUM_TIMESTEPS) data->globals.slkIndex0 = 0;
   if (data->globals.slkIndexLast>=SLK_NUM_TIMESTEPS) data->globals.slkIndexLast = 0;
#endif

   data->globals.time += data->globals.dt/2;

   data->reset();
#ifdef CORRECTED_GRADIENT
   data->neighboursGroup<calcInvM,ifSphOrSphBoundary>();
#endif
#ifdef VAR_H_CORRECTION 
   data->neighboursGroup<calcOmega,ifSphOrSphBoundary>();
#endif
   timeval t1,t2;
   gettimeofday(&t1,NULL);
#ifdef SLK
   data->neighboursGroup<calcGLeastSquares,ifSphOrSphBoundary>();
   data->syncParticlesBetweenProcs<vectTensor,readG,writeG>();
#endif
   data->neighbours<calcDddtDudt,ifSphOrSphBoundary>();
#ifdef LIQ_DEM
   data->syncParticlesBetweenProcs<vect,readF,writeF>();
   data->syncParticlesBetweenProcs<vect,readFb,writeFb>();
   data->syncParticlesBetweenProcs<vect,readFp,writeFp>();
   data->syncParticlesBetweenProcs<vect,readFv,writeFv>();
   data->neighboursGroupCoupling<calcFPI,ifDem>();
#ifndef FIXED_DEM
#ifndef NO_DEM_CONTACTS
   data->neighboursGroup<calcDEMContacts,ifDem>();
#ifndef LIQ_DEM_CUSTOM_WALL_CONTACTS
   data->neighboursGroup<calcDEMWallContacts,ifDem>();
#endif
#endif
#endif
#ifdef DEM_BLOCK
   vect avForce;
   avForce = 0.0;
   data->globals.tmp = 0;
   data->traverse<vect,calcAvForce,ifDem>(avForce);
   avForce = avForce/data->globals.tmp;
   data->traverse<vect,setAvForce,ifDem>(avForce);
#endif
   data->traverse<kick,ifDem>(0.5*data->globals.dt,0.5*data->globals.dt);
#endif
   gettimeofday(&t2,NULL);
   data->globals.wtTotalDddtCalc.tv_sec += t2.tv_sec-t1.tv_sec;
   data->globals.wtTotalDddtCalc.tv_usec += t2.tv_usec-t1.tv_usec;

   data->traverse<driftRest,ifSphOrSphBoundary>();


   
   data->traverse<calcEnergies,ifSphOrSphBoundary>();

}




