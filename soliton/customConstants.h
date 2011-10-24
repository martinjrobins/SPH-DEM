
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;

#include "parameters.h" 

const double RMIN[NDIM] = {0,0,0};
const double RMAX[NDIM] = {L1+L2+L3,1,H1+H2+H3};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {0,0,0};
const double DENS_DROP[NDIM] = {0,0,0};

const double DENS = 997.0479;
const double REFD = 997.0479;
//#define QUINTIC
#define WENDLAND
#define VISC_MORRIS
//#define HANN
//#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//#define CORRECTED_GRADIENT
//#define DENS_DIFFUSE
#define NO_ANTICLUMPING
#define CONST_H
//#define MORRIS_SPH_BOUNDARY
//#define CHECK_FOR_NAN
const double MAXTIME = 5.0;
const double DAMPTIME = 0.0;
const int OUTSTEP = 1000;
const int RESTART_EVERY_N_STEPS = 250;
const int REINIT_DENS_EVERY_N_STEPS = 500000000;

const double VISCOSITY = 8.92635148e-4/DENS;

const double VREF = sqrt(2.0*9.81*RMAX[2]);
const double VMAX = VREF;

const double REYNOLDS_NUMBER = VREF*(RMAX[0]-RMIN[0])/VISCOSITY;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;

const double MY_VAR_RES_DV = VMAX/8.0;

const int NY = 5;
const int NX = NY*RMAX[0];
const int NZ = NY*RMAX[2];
const double PSEP = (RMAX[1]-RMIN[1])/NY;
const int MAX_NUM_PARTICLES_PER_CPU = NX*NY*NZ;
const double BFAC = 1.00;
const double H = HFAC*PSEP;


const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


#define MIN_ALPHA 0.1
#endif
