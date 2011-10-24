
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS

const double PI = 3.14159265;

//next time I do this, RMIN and RMAX should correspond to the domain size and
//I should have a separate BMIN and BMAX for the box size
const double RMIN[NDIM] = {-1.0,-1.0};
const double RMAX[NDIM] = {1.0,1.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {0,0};
const int GHOST[2*NDIM] = {2,2,2,2};
const double DENS_DROP[NDIM] = {0.0,0.0};

const double DENS = 1000.0;
const double REFD = 1000.0;
#define WENDLAND
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
//#define DDENS_VARIANT
//#define BACKCOMPAT_READ
//#define MORRIS_SPH_BOUNDARY
//#define INCL_THERM_ENERGY

//#define PAIR_DBL
//#define FFT

const double OMEGA_E1 = 320;
const double OMEGA_E2 = -320;
const double CENTRE1[2] = {0,0.1};
const double CENTRE2[2] = {0,-0.1};
const double R0 = 0.1;

const int DAMPTIME = 0;
const double MAXTIME = 2.0;
const int OUTSTEP = 300;
const int RESTART_EVERY_N_STEPS = 100;
const int REINIT_DENS_EVERY_N_STEPS = 500000000;

const double VREF = 1.0;
const double REYNOLDS_NUMBER = 2500;
const double VMAX = 11.0;
const double VISCOSITY = VREF*0.5*(RMAX[0]-RMIN[0])/REYNOLDS_NUMBER;

const double CSFAC = 10.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;


//const double ALPHA = 0.1;
//const double H = Nmisc::viscToH(VISCOSITY,ALPHA,SPSOUND);
//const double PSEP = H/HFAC;
//const int NX = static_cast<int>((RMAX[0]-RMIN[0])/PSEP);
//const int NY = static_cast<int>((RMAX[1]-RMIN[1])/PSEP);

const int NX = 500;
const int NY = 500;
const int NZ = 0;
const double PSEP = (RMAX[0]-RMIN[0])/NX;
//const double H = HFAC*PSEP;
const double H = HFAC*(RMAX[0]-RMIN[0])/NX;
const double BFAC = 1.0;

const double APPROX_SMOOTH_GRID_SIZE = H/2;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);
const double MAX_NUM_PARTICLES_PER_CPU = NX*NY;

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
