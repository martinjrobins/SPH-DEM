
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;

const double RMIN[NDIM] = {0,0,0};
const double RMAX[NDIM] = {1.0,1.0,1.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {1,1,1};

const double DENS = 1000.0;
const double REFD = 1000.0;
const double PRESSURE_DROP[NDIM] = {0,0,0};
//#define QUINTIC
#define WENDLAND
//#define VISC_MORRIS
//#define HANN
#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//

#define LIQ_DEM
#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY

#define BACKGROUND_PRESSURE_FORCE
const double BGP_ACCEL[NDIM] = {0,0,9.81};


const double MAXTIME = 20;
const double DAMP_STEPS = 0;
const int OUTSTEP = 500;
const int RESTART_EVERY_N_STEPS = 10000;
const int REINIT_DENS_EVERY_N_STEPS = 2000000;

const double DEM_RADIUS = (RMAX[0]-RMIN[0])/80;
const double DEM_DENS = DENS*5;

const double REYNOLDS_NUMBER = 0.1;
const double VISCOSITY = sqrt((2.0/9.0)*(DEM_DENS-DENS)*9.81*pow(DEM_RADIUS,2)/(REYNOLDS_NUMBER*DENS));
const double VREF = REYNOLDS_NUMBER*VISCOSITY;
const double VMAX = VREF;

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;

//const double DENS_DROP[NDIM] = {pow(PRESSURE_DROP[0]/PRB+pow(DENS/REFD,7.0),1.0/7.0)*REFD-DENS,pow(PRESSURE_DROP[1]/PRB+pow(DENS/REFD,7.0),1.0/7.0)*REFD-DENS};
const double DENS_DROP[NDIM] = {0,0,0};

const double MY_VAR_RES_DV = VMAX/8.0;

const int NX = 20;
const int NY = 20;
const int NZ = 20;
const double PSEP = (RMAX[0]-RMIN[0])/NX;
//const double H = HFAC*PSEP;
const double H = HFAC*(RMAX[0]-RMIN[0])/NX;
const double BFAC = 1.00;

const double MAX_NUM_PARTICLES_PER_CPU = NX*NY*NZ*2.0;


const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
