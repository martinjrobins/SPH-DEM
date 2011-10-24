
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
const double BMAX[NDIM] = {2.0,1.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.3;

const int PERIODIC[NDIM] = {0,0};
const int GHOST[2*NDIM] = {0,0,0,0};
const double DENS_DROP[NDIM] = {0,0};



const double DENS = 1000.0;
const double REFD = 1000.0;

//#define RADIAL_BOUNDARY
#define WENDLAND
//#define QUINTIC
//#define HANN
//#define VISC_MORRIS
#define VISC_MONAGHAN
//#define VISC_CLEARY
//#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//#define SMOOTHING
//#define SPH_SMOOTHING
//const double EPSILON = 1.0;
//#define GRID_SMOOTHING
//#define SMOOTHED_VISC_VELOCITY
//#define SMOOTHED_VISC_VELOCITY_MLS
//#define SMOOTHED_VISC_VELOCITY_HAT

//#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY
//#define INCL_THERM_ENERGY
//#define SLK
//const int SLK_NUM_TIMESTEPS = 20;

//#define PAIR_DBL
//#define FFT

#define CORRECTED_GRADIENT


const int OUTSTEP = 500;
const int RESTART_EVERY_N_STEPS = 900;
const int REINIT_DENS_EVERY_N_STEPS = 1;
const double MAXTIME = 4.0;
const double DAMPTIME = 0.0;

const double VREF = 9.0;
const double REYNOLDS_NUMBER = 1000;
const double VMAX = VREF;
const double VISCOSITY = VREF*(BMAX[1]-BMIN[1])/REYNOLDS_NUMBER;

//#define MY_VAR_RES
//const double MY_VAR_RES_DV = VMAX/8.0;



const double CSFAC = 10.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;

const int NX = 20;
const double HR = 0.8;
const int NY = 80;
const int MAX_NUM_PARTICLES_PER_CPU = 4*NX*NY;
const double PSEP = (BMAX[1]-BMIN[1])*HR/NY;
const double BFAC = 1.0;
const double H = HFAC*PSEP;

const double APPROX_SMOOTH_GRID_SIZE = H/2;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

const double RMIN[NDIM] = {BMIN[0],BMIN[1]};
const double RMAX[NDIM] = {BMAX[0],BMAX[1]};

#define MIN_ALPHA 0.1
#endif
