
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS

const double PI = 3.14159265;

const double RMIN[NDIM] = {0,0};
const double RMAX[NDIM] = {20.0,20.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 7.0;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {1,0};

const double DENS = 1.0;
const double REFD = 0.95;
const double PRESSURE_DROP[NDIM] = {2.0,0};
//#define QUINTIC
#define WENDLAND
#define VISC_MORRIS
//#define HANN
//#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//

#define LIQ_DEM
#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY


const double MAXTIME = 20;
const double DAMP_STEPS = 50;
const int OUTSTEP = 500;
const int RESTART_EVERY_N_STEPS = 10000;
const int REINIT_DENS_EVERY_N_STEPS = 2000000;

const double VREF = 1.0;
const double VMAX = VREF;
const double VISCOSITY = 0.1/DENS;
const double REYNOLDS_NUMBER = VREF*(RMAX[0]-RMIN[0])/VISCOSITY;

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,GAMMA-1)*pow(SPSOUND,2)*REFD/GAMMA;

const double DENS_DROP[NDIM] = {pow(PRESSURE_DROP[0]/PRB+1.0,1.0/7.0)*REFD-REFD,pow(PRESSURE_DROP[1]/PRB+1.0,1.0/7.0)*REFD-REFD};

const double MY_VAR_RES_DV = VMAX/8.0;

const int NX = 30;
const int NY = 30;
const double PSEP = (RMAX[0]-RMIN[0])/NX;
//const double H = HFAC*PSEP;
const double H = HFAC*(RMAX[0]-RMIN[0])/NX;
const double BFAC = 1.00;

#define _2D_DEM
#define FIXED_DEM
#define ERGUN
const double DEM_RADIUS = 0.252313;
const double DEM_DENS = DENS*1;

const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
