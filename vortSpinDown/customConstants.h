
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS

const double PI = 3.14159265;

const double RMIN[NDIM] = {0,0};
const double RMAX[NDIM] = {1.0,1.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.95;

const int PERIODIC[NDIM] = {0,0};

const double DENS = 1000.0;
const double REFD = 990.0;
//#define QUINTIC
//#define VISC_MORRIS
#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME


const double MAXTIME = 13;
const int OUTSTEP = 100;
const int RESTART_EVERY_N_STEPS = 100;
const int REINIT_DENS_EVERY_N_STEPS = 200000000;

const double REYNOLDS_NUMBER = 420;
const double VREF = 1.0;
const double VMAX = sqrt(2.0)*(0.25/2.0);
const double VISCOSITY = VREF*(RMAX[0]-RMIN[0])/REYNOLDS_NUMBER;

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
//const double H = HFAC*PSEP;
const double H = HFAC*(RMAX[0]-RMIN[0])/NX;


const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
