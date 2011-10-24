
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include "misc.h"
#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS


const double RMIN[NDIM] = {-1.0,-1.0};
const double RMAX[NDIM] = {1.0,1.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

#define GAMMA 1.4
#define HFAC  1.95
#define PI    3.14159265

const double PMIN[NDIM] = {-100,-100};
const double PMAX[NDIM] = {100,100};

#define DENS 1000
#define REFD 1000
//#define QUINTIC
#define VISC_MORRIS
//#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
#define REFERENCE_FRAME

//#define PAIR_DBL
//#define FFT
const double PAIR_DBL_R_MULT2 = pow(1.2,2);
const double PAIR_DBL_MAX_TIME = 100;
const double PAIR_DBL_BINS = 10;
const int PAIR_DBL_SAMPLES_PER_OUTSTEP = 1000;

const double MAXTIME = 1;
const int OUTSTEP = 100;
const int RESTART_EVERY_N_STEPS = 50;
const int REINIT_DENS_EVERY_N_STEPS = 200000000;


const double VREF = 1;
const double REYNOLDS_NUMBER = 5000;
const double VISCOSITY = VREF*0.5*(RMAX[0]-RMIN[0])/REYNOLDS_NUMBER;

//const double VMAX = 3.00;
const double VMAX = VREF;
const double CSFAC = 10.0;
const double SPSOUND = CSFAC*VMAX;
//const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;
const double PRB = pow(SPSOUND,2)*REFD/7.0;


//const double ALPHA = 0.1;
//const double H = Nmisc::viscToH(VISCOSITY,ALPHA,SPSOUND);
//const double PSEP = H/HFAC;
//const int NX = static_cast<int>((RMAX[0]-RMIN[0])/PSEP);
//const int NY = static_cast<int>((RMAX[1]-RMIN[1])/PSEP);

const int NX = 20;
const int NY = 20;
const double PSEP = (RMAX[0]-RMIN[0])/NX;
//const double H = HFAC*PSEP;
const double H = HFAC*(RMAX[0]-RMIN[0])/NX;


const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
