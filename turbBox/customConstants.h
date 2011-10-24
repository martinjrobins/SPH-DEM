
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS

//#define SPIN
#define FORCING
const double PI = 3.14159265;

const double FORCE_AMP = 6.0;
const double FORCE_TAU = 0.01683;
//const double FORCE_TAU = 1.0;

//next time I do this, RMIN and RMAX should correspond to the domain size and
//I should have a separate BMIN and BMAX for the box size
const double RMIN[NDIM] = {-1.0,-1.0};
const double RMAX[NDIM] = {1.0,1.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.95;

const int PERIODIC[NDIM] = {0,0};
const int GHOST[2*NDIM] = {0,0,0,0};
const double DENS_DROP[NDIM] = {0.0,0.0};

const double DENS = 1000.0;
const double REFD = 1000.0;
#define WENDLAND
//#define QUINTIC
//#define HANN
//#define VISC_MORRIS
#define VISC_MONAGHAN
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
#define DDENS_VARIANT
//#define BACKCOMPAT_READ
//#define MORRIS_SPH_BOUNDARY
//#define INCL_THERM_ENERGY

//#define PAIR_DBL
//#define FFT
#ifdef SPIN
const double BOX_GLOB_ANGVEL = 1;
const double BOX_PERT_ANGVEL = 0.75;
const double BOX_PERT_FREQ = 1.0/(2.0*PI);
#else
const double BOX_GLOB_ANGVEL = 0;
const double BOX_PERT_ANGVEL = 0.0;
const double BOX_PERT_FREQ = 1.0/(2.0*PI);
#endif
#ifdef SPIN
const double BOX_GLOB_TIME = 0.1/BOX_PERT_FREQ;
#else 
const double BOX_GLOB_TIME = 0.0/BOX_PERT_FREQ;
#endif

const double NOREVS = 10;
const int DAMP_STEPS = 0;
const double MAXTIME = BOX_GLOB_TIME+ NOREVS/BOX_PERT_FREQ;
const int OUTSTEP = 200;
const int RESTART_EVERY_N_STEPS = 500;
const int REINIT_DENS_EVERY_N_STEPS = 500000000;

#ifdef SPIN
const double VREF = BOX_PERT_ANGVEL*sqrt(pow(0.5*(RMAX[0]-RMIN[0]),2)+pow(0.5*(RMAX[1]-RMIN[1]),2));
const double MAX_VRAND = 1.0;
const double REYNOLDS_NUMBER = 5000;
const double VMAX = 2.60;
const double VISCOSITY = VREF*0.5*(RMAX[0]-RMIN[0])/REYNOLDS_NUMBER;
#else
#ifdef FORCING
const double VREF = 0.5;
const double KE_REF = 1.0;
//const double REYNOLDS_NUMBER = VREF*0.5*(RMAX[0]-RMIN[0])/0.0005;
const double VMAX = 2.0;
const double VISCOSITY = 0.0005;
//const double VISCOSITY = 0.000001;
//const double VISCOSITY = 0.00000000;
const double REYNOLDS_NUMBER = VREF*0.5*(RMAX[0]-RMIN[0])/VISCOSITY;
#else
const double VREF = 1.0;
const double KE_REF = 0.5*pow(VREF,2)*(RMAX[1]-RMIN[1])*(RMAX[0]-RMIN[0])*DENS;
const double REYNOLDS_NUMBER = 1500;
const double VMAX = 3.00;
const double VISCOSITY = VREF*0.5*(RMAX[0]-RMIN[0])/REYNOLDS_NUMBER;
#endif
#endif

const double CSFAC = 10.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;


//const double ALPHA = 0.1;
//const double H = Nmisc::viscToH(VISCOSITY,ALPHA,SPSOUND);
//const double PSEP = H/HFAC;
//const int NX = static_cast<int>((RMAX[0]-RMIN[0])/PSEP);
//const int NY = static_cast<int>((RMAX[1]-RMIN[1])/PSEP);

const int NX = 100;
const int NY = 100;
const int NZ = 0;
const double PSEP = (RMAX[0]-RMIN[0])/NX;
//const double H = HFAC*PSEP;
const double H = HFAC*(RMAX[0]-RMIN[0])/NX;
const double BFAC = 1.0;

const double APPROX_SMOOTH_GRID_SIZE = H/2;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
