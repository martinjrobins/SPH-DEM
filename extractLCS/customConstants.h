
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;
const double DENS = 1;
const double REFD = DENS;

const double RMIN[NDIM] = {-0.082,-0.082,-0.01};
const double RMAX[NDIM] = {0.082,0.082,0.155};
const double BMIN[NDIM] = {-0.082,-0.082,-0.01};
const double BMAX[NDIM] = {0.082,0.082,0.155};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 7.0;
const double HFAC = 1.5;

//#define GHOST_BOUNDARY
const int PERIODIC[NDIM] = {0,0,0};
const int GHOST[2*NDIM] = {0,0,0,0,0,0};
const double DENS_DROP[NDIM] = {0,0,0};

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
//#define VAR_H_CORRECTION

//#define SETTLE
//#define SPHBOUNDARY
//#define LIQ_DEM
//#define LIQ_DEM_TEST
//#define LIQ_DEM_SEPARATE_DRAGS
//#define LIQ_DEM_SIMPLE_DRAG
//#define LIQ_DEM_ADDED_MASS
//#define LIQ_DEM_ONE_WAY_COUPLE
//#define LIQ_DEM_CUSTOM_WALL_CONTACTS
//#define LIQ_DEM_DDDT_VER2
//#define LINEAR
//#define MANY_PARTICLES
//#define DENS_DIFFUSE

#define NO_ANTICLUMPING
#define CONST_H
//#define MORRIS_SPH_BOUNDARY

//const double Ar = 8*pow(DEM_RADIUS,3)*(DEM_DENS-DENS)*9.81/VISCOSITY

const double VREF = 1.0;
const double VMAX = VREF;
const double VISCOSITY = 1.0;
#define DEFORMATION_MATRIX
#define HALLER_LCS
const double REYNOLDS_NUMBER = VREF*0.5*(RMAX[0]-RMIN[0])/VISCOSITY;

const double MAXTIME = 2*PI;
const int DAMPTIME = 0.0;
const int OUTSTEP = 2;
const int RESTART_EVERY_N_STEPS = 250;
const int REINIT_DENS_EVERY_N_STEPS = 50000000;
const int REINIT_DENS_AT_N_STEPS = 50000000;
//#define REINIT_DENS_MLS



//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,GAMMA-1.0)*pow(SPSOUND,2)*REFD/GAMMA;


const double MY_VAR_RES_DV = VMAX/8.0;

const int NX = 36*2;
const int NY = NX;
//const double PSEP = (BMAX[0]-BMIN[0])/(NX);
const double PSEP = 0.002;
const int NZ = 71;
const int MAX_NUM_PARTICLES_PER_CPU = NX*NY*NZ;
const double BFAC = 1.00;
const double H = HFAC*PSEP;


const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


#define MIN_ALPHA 0.1
#endif
