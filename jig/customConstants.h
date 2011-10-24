
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
const double HFAC = 1.95;

const int PERIODIC[NDIM] = {1,1};

const double DENS = 1000.0;
const double REFD = 995.0;
//#define QUINTIC
//#define PRICE2
//const double PRICE2_ALPHA = 1.5;
//const double PRICE2_BETA = 0.7;
//#define MY_CUBIC
//const double MY_CUBIC_ALPHA = 1.0/HFAC;
#define WENDLAND
//#define VISC_MORRIS
#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//#define SMOOTHING
//#define SPH_SMOOTHING
//const double EPSILON = 1.5;
//#define GRID_SMOOTHING

#define NO_ANTICLUMPING
#define CONST_H
//#define MORRIS_SPH_BOUNDARY
//#define INCL_THERM_ENERGY

//#define PAIR_DBL
//#define FFT

const int OUTSTEP = 1000;
const int RESTART_EVERY_N_STEPS = 200;
const int REINIT_DENS_EVERY_N_STEPS = 500000000;

const double VREF = 1.0;
const double REYNOLDS_NUMBER = 1500;
const double VMAX = 1.0;
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

const int NX = 50;
const int NY = 50;
const double PSEP = (RMAX[0]-RMIN[0])/NX;
const double MAXTIME = 1000*PSEP/(2*VMAX);
const double JIG_A = 0.5*PSEP;
const double JIG_FREQ = 2*VMAX/PSEP;
const double H = HFAC*PSEP;

const double APPROX_SMOOTH_GRID_SIZE = 2*H;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
