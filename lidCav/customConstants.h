
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS

const double BMIN[NDIM] = {-1.0,-1.0};
const double BMAX[NDIM] = {1.0,1.0};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

#define GAMMA 1.4
#define HFAC  1.95
#define PI    3.14159265

const int PERIODIC[NDIM] = {1,0};

#define DENS 1000.0
#define REFD 1000.0
//#define QUINTIC
//#define VISC_MORRIS
#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME

const double VREF = 1.0;
const double NOREVS = 40.0;
const double MAXTIME = NOREVS*(BMAX[0]-BMIN[0])/VREF;
const int OUTSTEP = 2000;
const int RESTART_EVERY_N_STEPS = 100;
const int REINIT_DENS_EVERY_N_STEPS = 200000000;

const double REYNOLDS_NUMBER = 1000;
const double VISCOSITY = VREF*(BMAX[0]-BMIN[0])/REYNOLDS_NUMBER;

const double VMAX = VREF;
const double CSFAC = 10.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;

#define MY_VAR_RES
const double MY_VAR_RES_DV = VMAX/2.0;

//const double ALPHA = 0.1;
//const double H = Nmisc::viscToH(VISCOSITY,ALPHA,SPSOUND);
//const double PSEP = H/HFAC;
//const int NX = static_cast<int>((RMAX[0]-RMIN[0])/PSEP);
//const int NY = static_cast<int>((RMAX[1]-RMIN[1])/PSEP);

const int NX = 50;
const int NY = 50;
const double PSEP = (BMAX[0]-BMIN[0])/NX;
//const double H = HFAC*PSEP;
const double H = HFAC*(BMAX[0]-BMIN[0])/NX;
const double GRIDSEP = PSEP;
#define MIN_ALPHA 0.1
const double RMIN[NDIM] = {BMIN[0]-3.5*PSEP,BMIN[1]-3.5*PSEP};
const double RMAX[NDIM] = {BMAX[0]+3.5*PSEP,BMAX[1]+3.5*PSEP};


const double APPROX_SMOOTH_GRID_SIZE = H/2;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);
#endif
