
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;

const double RMIN[NDIM] = {0,0,0};
const double RMAX[NDIM] = {1,1,1};
const double BMIN[NDIM] = {0,0,0};
const double BMAX[NDIM] = {1,1,1};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {1,1,1};
const double DENS_DROP[NDIM] = {0,0,0};

const double DENS = 1000.0;
const double REFD = 1000.0;
//#define QUINTIC
#define WENDLAND
#define VISC_MORRIS
//#define HANN
//#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//

//#define LIQ_DEM
#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY
#define FIXED_DEM

//#define CHECK_FOR_NAN
const double POROSITY = 0.65;
const double MAXTIME = 1.0;
const int DAMP_STEPS = 0;
const int OUTSTEP = 100;
const int RESTART_EVERY_N_STEPS = 50;
const int REINIT_DENS_EVERY_N_STEPS = 2000000;

const double VREF = 1.0;
const double VMAX = VREF;
const double REYNOLDS_NUMBER = 1.0;
const double VISCOSITY = VREF*(BMAX[0]-BMIN[0])/REYNOLDS_NUMBER;

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;


const double MY_VAR_RES_DV = VMAX/8.0;

const int NX = 10;
const int NY = NX;
const int NZ = NX;
const double PSEP = (BMAX[0]-BMIN[0])/(NX);
const int MAX_NUM_PARTICLES_PER_CPU = 4.0*NX*NY*NZ;
const double BFAC = 1.00;
//const double H = HFAC*PSEP;
const double H = HFAC*PSEP;

const double DEM_RADIUS = (RMAX[0]-RMIN[0])/25;
const double DEM_DENS = DENS*5.0;

const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
