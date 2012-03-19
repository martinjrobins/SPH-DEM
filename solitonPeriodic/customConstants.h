
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;

#include "parameters.h" 

const int STARTTAG = 1;
const int ENDTAG = 2;
const int NY = (PY*2+1)*(RY*2+1);
const double PSEP = WIDTH/NY;
const int NX = (L1+L2+L3+L4)/PSEP;
const int NZ = int((H1+H2+H3)/PSEP);
const int NXP = NX-L4/PSEP;


const double LL = L1+L2+L3;
const double TL = LL + (L4+L5)*(RY*2+1);
const int ENDPERIODICCPU = int((LL/TL)*NCPU+0.5)-1;

const double RMIN[NDIM] = {0,0,0};
const double RMAX[NDIM] = {L1+L2+L3+L4+L5,WIDTH,WALLH*1.5};


//#define MODIFY_BFORCE_WITH_STILL_LEVEL
//const double STILL_LEVEL = H1;

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {0,1,0};
const int GHOST[2*NDIM] = {1,0,0,0,1,0};
const double DENS_DROP[NDIM] = {0,0,0};

const double DENS = 997.0479;
const double REFD = 997.0479;
//#define QUINTIC
#define WENDLAND
//#define VISC_MORRIS
//#define HANN
//#define VISC_MONAGHAN
#define VISC_ARTIFICIAL
const double ALPHA_ARTIFICIAL = 0.1;
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//#define CORRECTED_GRADIENT
//#define DENS_DIFFUSE
#define NO_ANTICLUMPING
#define CONST_H
#define SLIP_BOUNDARIES
//#define CORRECTED_GRADIENT
//#define MORRIS_SPH_BOUNDARY
//#define CHECK_FOR_NAN
const double DAMPTIME = WALLUP/2.0;
const int RESTART_EVERY_N_STEPS = 100;
const int REINIT_DENS_EVERY_N_STEPS = 20000000;


const double VREF = sqrt(2.0*9.81*RMAX[2]);
const double VMAX = VREF;

//const double VISCOSITY = VREF*(RMAX[2]-RMIN[2])/REYNOLDS_NUMBER;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;

const double MY_VAR_RES_DV = VMAX/8.0;

const int MAX_NUM_PARTICLES_PER_CPU = ((L1+L2+L3)/(PSEP*(ENDPERIODICCPU+1)))*(RY*2+1)*3*NZ*1.2;
const double BFAC = 1.00;
const double H = HFAC*PSEP;


const double VISCOSITY = (15.0/112.0)*0.1*H*SPSOUND;
const double REYNOLDS_NUMBER = VREF*(RMAX[2]-RMIN[2])/VISCOSITY;

const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


#define MIN_ALPHA 0.1
#endif
