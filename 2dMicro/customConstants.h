
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS

const double PI = 3.14159265;


//next time I do this, RMIN and RMAX should correspond to the domain size and
//I should have a separate BMIN and BMAX for the box size
const double BMIN[NDIM] = {0.0,0.0};
const double BMAX[NDIM] = {2.857,1.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.3;

const int PERIODIC[NDIM] = {1,0};

const double DENS = 1000.0;
const double REFD = 980.0;

//#define RADIAL_BOUNDARY
//#define WENDLAND
//#define QUINTIC
//#define HANN
#define VISC_MORRIS
//#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//#define SMOOTHING
//#define SPH_SMOOTHING
//const double EPSILON = 1.0;
//#define GRID_SMOOTHING
//#define SMOOTHED_VISC_VELOCITY
//#define SMOOTHED_VISC_VELOCITY_MLS
//#define SMOOTHED_VISC_VELOCITY_HAT

//#define NO_ANTICLUMPING
#define CONST_H
//#define MORRIS_SPH_BOUNDARY
//#define INCL_THERM_ENERGY

//#define PAIR_DBL
//#define FFT

//this contains H,R and WRATIO
#include "mixerParameters.h"
#define HALF_COURANT

const double THETA = 22.5*PI/180.0;
const double H1 = 0.2*(BMAX[1]-BMIN[1]);
//#define NO_TOP
const double H2 = HRATIO * H1;
const double WW1 = 0.5*(BMAX[0]-BMIN[0]);
const double WW2 = WRATIO * WW1;
const double V1 = 1.0;
const double V2 = VRATIO*V1;
const double NOPERIODS = 19.0;
const int OUTSTEP = 20;
const int RESTART_EVERY_N_STEPS = 900;
const int REINIT_DENS_EVERY_N_STEPS = 500000000;

const double VREF = V1;
const double REYNOLDS_NUMBER = 10;
const double VMAX = V1;
const double VISCOSITY = VREF*(BMAX[1]-BMIN[1])/REYNOLDS_NUMBER;

const double CSFAC = 10.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;

const int NX = 100;
const double PSEP = (BMAX[0]-BMIN[0])/NX;
const int NY = int((BMAX[1]-BMIN[1])/PSEP);
const double H = HFAC*PSEP;

const double APPROX_SMOOTH_GRID_SIZE = H/2;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


const double RMIN[NDIM] = {0.0-1.1*H,0.0};
const double RMAX[NDIM] = {2.857+1.1*H,1.0};

const double DAMP = 0.5*(RMAX[0]-RMIN[0])/V1;
const double SPEEDUP = 0.5*(RMAX[0]-RMIN[0])/V1;
const double MAXTIME = NOPERIODS*(RMAX[0]-RMIN[0])/V1 + DAMP + SPEEDUP;
#define MIN_ALPHA 0.1
#endif
