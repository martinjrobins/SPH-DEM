
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS

const double PI = 3.14159265;

const double RMIN[NDIM] = {0,0};
const double RMAX[NDIM] = {1,1};
const double BMIN[NDIM] = {0,0};
const double BMAX[NDIM] = {1,1};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 7.0;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {1,1};
const int GHOST[2*NDIM] = {0,0,0,0};

const double DENS = 1000.0;
const double REFD = 1000.0;
//#define QUINTIC
#define WENDLAND
#define VISC_MORRIS
//#define VISC_ARTIFICIAL
//const double ALPHA_ARTIFICIAL = 0.1;
//#define HANN
//#define VISC_MONAGHAN
//#define VISC_CLEARY
//#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//#define CORRECTED_GRADIENT
#define VAR_H_CORRECTION
//#define CONST_H

//#define SETTLE
//#define SPHBOUNDARY
#define LIQ_DEM
//#define LIQ_DEM_SEPARATE_DRAGS
//#define LIQ_DEM_DDDT_VER2
#define _2D_DEM
//#define LIQ_DEM_SIMPLE_DRAG
//#define LIQ_DEM_ADDED_MASS
//#define LIQ_DEM_ONE_WAY_COUPLE
#define LINEAR
//#define DENS_DIFFUSE

#include "parameters.h" 

const double DEM_RADIUS = (RMAX[0]-RMIN[0])/80;
const double DEM_K = 2.5*10e5;
#ifdef _2D_DEM
const double DEM_VOL = PI*pow(DEM_RADIUS,2);
#else
const double DEM_VOL = (4.0/3.0)*PI*pow(DEM_RADIUS,3);
#endif
const double DEM_MASS = DEM_VOL*DEM_DENS;
const double DEM_GAMMA = 2.0;
const double DEM_MIN_REDUCED_MASS = 0.5*DEM_MASS;
#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY
#define FIXED_DEM
//#define LIQ_DEM_DENS_NORMAL

//#define CHECK_FOR_NAN
//const double POROSITY = 0.90;
const double MAXTIME = 10.0;
const double TIME_DROP_PARTICLE = 0.0;
const int DAMPTIME = 0;
const int OUTSTEP = 300;
const int RESTART_EVERY_N_STEPS = 250;
const int REINIT_DENS_EVERY_N_STEPS = 50;
const int REINIT_DENS_AT_N_STEPS = 50000000;
//#define REINIT_DENS_MLS


//const double Ar = 8*pow(DEM_RADIUS,3)*(DEM_DENS-DENS)*9.81/VISCOSITY

const double REYNOLDS_NUMBER = 0.01;
const double Ar = REYNOLDS_NUMBER*18;
const double VISCOSITY = sqrt(8*pow(DEM_RADIUS,3)*DENS*(DEM_DENS-DENS)*9.81/Ar)/DENS;
const double VREF = 0.5*VISCOSITY*REYNOLDS_NUMBER/DEM_RADIUS;
const double VMAX = 20.0;

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(SPSOUND,2)*DENS/(pow(DENS/REFD,GAMMA)*(GAMMA+1)-1);

const double DENS_DROP[NDIM] = {pow(10*REFD/PRB+1,1.0/GAMMA)*REFD-REFD,0};

const double MY_VAR_RES_DV = VMAX/8.0;

const int NY = NX;
const double PSEP = (BMAX[0]-BMIN[0])/(NX);
const int MAX_NUM_PARTICLES_PER_CPU = 2.0*NX*NY;
const int MAX_NUM_DEM_PARTICLES = 2.0*NX*NY;
const double BFAC = 1.00;
//const double H = HFAC*PSEP;
const double H = HFAC*PSEP;

const double LIQ_DEM_COUPLING_RADIUS = 1.9999*H/pow(POROSITY,1.0/NDIM);

const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


#define MIN_ALPHA 0.1
#endif
