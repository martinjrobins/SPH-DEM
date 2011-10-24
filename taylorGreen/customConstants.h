
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
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {1,1};

const double DENS = 1000.0;
const double REFD = 1000.0;
//#define QUINTIC
#define WENDLAND
//#define VISC_MORRIS
//#define HANN
#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//

//#define SMOOTHING
const double EPSILON = 0.9;
//#define DIRECT_SMOOTHING
//#define MY_VAR_RES
//#define SLK
//const int SLK_NUM_TIMESTEPS = 20;

//#define DS_ITERATE
//#define DENS_DIFFUSE
//const double SMOOTH_HFAC_MULT = 1.0;
//#define SPH_SMOOTHING
//const double EPSILON = 1.5;
//#define GRID_SMOOTHING

#define LIQ_DEM
const double DEM_RADIUS = 1.0;
const double DEM_DENS = 1000;

#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY
//#define INCL_THERM_ENERGY



const double MAXTIME = 6;
const int OUTSTEP = 1000;
const int RESTART_EVERY_N_STEPS = 100;
const int REINIT_DENS_EVERY_N_STEPS = 2000000;
const int DAMP_STEPS = 0;

const double NVORT = 2;
const double REYNOLDS_NUMBER = 100;
const double VREF = 1.0;
const double VMAX = VREF;
const double VISCOSITY = VREF*(RMAX[0]-RMIN[0])/REYNOLDS_NUMBER;

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;


const double MY_VAR_RES_DV = VMAX/8.0;

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
const double BFAC = 1.0;

const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
