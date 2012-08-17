#include "sphIncompress.h"

void CsphIncompress::start() {
   data->traverse<driftRAndKick,ifSphOrSphBoundaryOrImmersedDem>();
#ifdef LIQ_DEM
   data->traverse<driftR,ifDem>();
#endif


#ifdef SLK
   data->globals.slkIndex0++;
   data->globals.slkIndexLast++;
   if (data->globals.slkIndex0>=SLK_NUM_TIMESTEPS) data->globals.slkIndex0 = 0;
   if (data->globals.slkIndexLast>=SLK_NUM_TIMESTEPS) data->globals.slkIndexLast = 0;
#endif
   data->traverse<driftRest,ifSphOrSphBoundaryOrImmersedDem>();
   data->globals.time += data->globals.dt/2;

   //if ((data->globals.sphStep < DAMP_STEPS)||(data->globals.sphStep%REINIT_DENS_EVERY_N_STEPS==0)) {
   if ((data->globals.sphStep%REINIT_DENS_EVERY_N_STEPS==0)||(data->globals.sphStep == REINIT_DENS_AT_N_STEPS)) {
      data->reset();
#ifdef MY_VAR_RES
      data->neighboursGroup<calcGradVLeastSquares,ifSphOrSphBoundaryOrImmersedDem>();
      data->syncParticlesBetweenProcs<vect,readGradV0,writeGradV0>();
      data->syncParticlesBetweenProcs<vect,readGradV1,writeGradV1>();
      //data->syncParticlesBetweenProcs<vect,readGradV2,writeGradV2>();
#endif
      //if (data->globals.mpiRank==0) cout<<"\tRe-initialising density...sphStep="<<data->globals.sphStep<<endl;
      for (int i=0;i<10;i++) {
         data->traverse<calcHFromDens,ifSphOrSphBoundaryOrImmersedDemOrGhost>();
         data->syncParticlesBetweenProcs<double,readH,writeH>();
         data->neighboursGroup<calcDensity,ifSphOrSphBoundaryOrImmersedDem>();
         data->syncParticlesBetweenProcs<double,readDens,writeDens>();
      }
   }

   data->traverse<initSumsMiddle,ifSphOrSphBoundaryOrImmersedDem>();
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

   timeval t1,t2;
   gettimeofday(&t1,NULL);
#ifdef LIQ_DEM
   data->neighboursGroupCoupling<calcPorosity,ifSphOrSphBoundaryOrImmersedDem>();
#ifndef LIQ_DEM_ONE_WAY_COUPLE
#ifdef LIQ_DEM_SEPARATE_DRAGS
   data->neighboursGroupCoupling<calcSPHDrag,ifSphOrSphBoundaryOrImmersedDem>();
#else
   data->traverse<initShepSum,ifDemOrGhost>();
   data->neighboursGroupCoupling<calcShepSum,ifSphOrSphBoundaryOrImmersedDem>();
   data->reverseSyncParticlesBetweenProcs<double,readShepSum,sumShepSum>();
   data->syncParticlesBetweenProcs<double,readShepSum,writeShepSum>();
   data->neighboursGroupCoupling<calcDEMtoSPHDrag,ifSphOrSphBoundaryOrImmersedDem>();
#endif
#endif
   data->traverse<finalisePorosity,ifSphOrSphBoundaryOrImmersedDem>();
   data->traverse<setPorosityForBoundary,ifBoundary>();
   data->syncParticlesBetweenProcs<double,readPorosity,writePorosity>();

#ifdef VAR_H_CORRECTION2
   /*
   data->traverse<initGradH,ifSphOrSphBoundaryOrDemOrGhost>();
   data->neighboursGroupCoupling<calcGradPorosity,ifDem>();
   data->reverseSyncParticlesBetweenProcs<vect,readGradH,sumGradH>();
   */
   data->neighboursGroupCoupling<calcGradPorosity2,ifSphOrSphBoundaryOrImmersedDem>();
#endif
#endif

#ifdef CORRECTED_GRADIENT
   data->neighboursGroup<calcInvM,ifSphOrSphBoundaryOrImmersedDem>();
#endif
#ifdef VAR_H_CORRECTION 
   data->neighboursGroup<calcOmega,ifSphOrSphBoundaryOrImmersedDem>();
   data->syncParticlesBetweenProcs<double,readOmega,writeOmega>();
#endif
#ifdef VAR_H_CORRECTION2
#ifdef LIQ_DEM
   data->traverse<finaliseGradH,ifSphOrSphBoundaryOrImmersedDem>();
#else
   data->neighboursGroup<calcGradH,ifSphOrSphBoundaryOrImmersedDem>();
#endif
   data->syncParticlesBetweenProcs<vect,readGradH,writeGradH>();
#endif

   data->traverse<calcPressSpsoundPdr2,ifSphOrSphBoundaryOrImmersedDemOrGhost>();	

#ifdef SLK
   data->neighboursGroup<calcGLeastSquares,ifSphOrSphBoundaryOrImmersedDem>();
   data->syncParticlesBetweenProcs<vectTensor,readG,writeG>();
#endif
#ifdef SMOOTHED_VISC_VELOCITY
   data->neighboursGroupScaled<smoothViscV,ifSphOrSphBoundaryOrImmersedDem>();
   data->syncParticlesBetweenProcs<vect,readViscV,writeViscV>();
#endif
   //cout <<"rank = "<<data->globals.mpiRank<<"before force "<<endl;
   data->neighbours<calcForce,ifSphOrSphBoundary>();
   //cout <<"rank = "<<data->globals.mpiRank<<"after force"<<endl;
#ifdef BACKGROUND_PRESSURE_FORCE
   data->traverse<addBackgroundPressure,ifSphOrSphBoundaryOrImmersedDem>();
#endif
#ifdef LIQ_DEM
#ifdef LIQ_DEM_TEST
   data->traverse<addGradPorosityTerm,ifSphOrSphBoundary>();
#endif
#ifdef LIQ_DEM_VARIABLE_TIMESTEP
   data->globals.newDt = min(data->globals.newDt,data->globals.demDt);
#else
   data->globals.newDt = min(data->globals.newDt,liqDemCondition());
#endif
#ifndef FIXED_DEM
#ifndef NO_DEM_CONTACTS
   if (data->globals.time > TIME_DROP_PARTICLE) {
      data->globals.newDt = min(data->globals.newDt,demCondition());
   }
#endif
#endif
#endif
   data->traverse<accelCondition,ifSphOrSphBoundaryOrImmersedDemOrDem>(); 
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
   

   data->traverse<initSumsEnd,ifSphOrSphBoundaryOrImmersedDem>();
#ifdef LIQ_DEM
   data->traverse<driftRAndKick,ifDem>();
   data->traverse<initSumsMiddle,ifDem>();
#endif
}

void CsphIncompress::end() {

#ifdef LIQ_DEM
   data->globals.demDt = data->globals.maxdt;
#endif

#ifndef SMOOTHING
   data->traverse<driftR,ifSphOrSphBoundaryOrImmersedDem>();
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
   data->neighboursGroup<calcInvM,ifSphOrSphBoundaryOrImmersedDem>();
#endif
#ifdef VAR_H_CORRECTION 
   data->neighboursGroup<calcOmega,ifSphOrSphBoundaryOrImmersedDem>();
#endif
#ifdef VAR_H_CORRECTION2
#ifdef LIQ_DEM
   data->neighboursGroupCoupling<calcGradPorosity2,ifSphOrSphBoundaryOrImmersedDem>();
   data->traverse<finaliseGradH,ifSphOrSphBoundaryOrImmersedDem>();
#else
   data->neighboursGroup<calcGradH,ifSphOrSphBoundaryOrImmersedDem>();
#endif
   data->syncParticlesBetweenProcs<vect,readGradH,writeGradH>();
#endif
#ifdef LIQ_DEM_DDDT_VER2
   data->neighboursGroupCoupling<calcPorosity2,ifSphOrSphBoundaryOrImmersedDem>();
   data->traverse<finalisePorosity,ifSphOrSphBoundaryOrImmersedDem>();
   data->syncParticlesBetweenProcs<double,readPorosity,writePorosity>();
   data->syncParticlesBetweenProcs<double,read_dPorositydt,write_dPorositydt>();
#endif


   timeval t1,t2;
   gettimeofday(&t1,NULL);
#ifdef SLK
   data->neighboursGroup<calcGLeastSquares,ifSphOrSphBoundaryOrImmersedDem>();
   data->syncParticlesBetweenProcs<vectTensor,readG,writeG>();
#endif
   data->neighbours<calcDddtDudt,ifSphOrSphBoundaryOrImmersedDem>();
#ifdef LIQ_DEM
   //data->traverse<correctDddt,ifSphOrSphBoundary>();
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

   data->traverse<driftRest,ifSphOrSphBoundaryOrImmersedDem>();


   
   data->traverse<calcEnergies,ifSphOrSphBoundaryOrImmersedDem>();

}




