
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS
#define MAXTIME 1.0


const double RMIN[NDIM] = {0.0,-0.00010};
const double RMAX[NDIM] = {0.0025,0.00110};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

#define GAMMA 1.4
#define HFAC  1.95
#define PI    3.14159265

const int NX = 50;
const int NY = 24;

const int MAX_NUM_PARTICLES_PER_CPU = 4*NX*NY;

const int PERIODIC[NDIM] = {1,0};
const double DENS_DROP[NDIM] = {0,0};

#define DENS 1000.0
#define REFD 1000.0
#define VREF 0.0025
#define CSFAC 10

//#define QUINTIC
//#define VISC_MORRIS
#define VISC_MONAGHAN
//#define VISC_CLEARY
//#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//#define CORRECTED_GRADIENT

const double PSEP = (RMAX[1]-RMIN[1])/NY;
const double BFAC = 1.0;
const double H = HFAC*PSEP;
const double PRB = pow(CSFAC*VREF,2)*REFD/7.0;
const double SPSOUND = sqrt(7.0*PRB/REFD);
const double GRIDSEP = PSEP;



const int DAMPTIME = 0.0;
const int OUTSTEP = 100;
const int RESTART_EVERY_N_STEPS = 1000;
const int REINIT_DENS_EVERY_N_STEPS = 20000;

const double REYNOLDS_NUMBER = 1.0125;
const double VISCOSITY = VREF*0.001/REYNOLDS_NUMBER;
//#define MY_VAR_RES
//const double MY_VAR_RES_DV = VREF/2.0;
//const double VISCOSITY = 0.0000125;
//const double REYNOLDS_NUMBER = VREF*0.001/VISCOSITY;

const double APPROX_SMOOTH_GRID_SIZE = H/2;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);


#define MIN_ALPHA 0.1
#endif
