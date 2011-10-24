
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;

const double RMIN[NDIM] = {0,0,0};
const double RMAX[NDIM] = {1,1,1.5};
const double BMIN[NDIM] = {0,0,0};
const double BMAX[NDIM] = {1,1,1.5};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {1,1,0};
const int GHOST[2*NDIM] = {0,0,0,0,1,0};
const double DENS_DROP[NDIM] = {0,0,0};

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
//#define CORRECTED_GRADIENT

//#define SETTLE
//#define SPHBOUNDARY
#define LIQ_DEM
//#define LIQ_DEM_SIMPLE_DRAG
//#define LIQ_DEM_ADDED_MASS
//#define LIQ_DEM_ONE_WAY_COUPLE
#define LINEAR
//#define NO_DEM_CONTACTS
//#define MANY_PARTICLES
//#define DENS_DIFFUSE

const double POROSITY = 0.6;
const double DEM_DENS = DENS*2.0;
const double DEM_RADIUS = (RMAX[0]-RMIN[0])/200;
const double DEM_K = 2.5*10e5;
const double DEM_VOL = (4.0/3.0)*PI*pow(DEM_RADIUS,3);
const double DEM_MASS = DEM_VOL*DEM_DENS;
const double DEM_GAMMA = 2.0;
const double DEM_MIN_REDUCED_MASS = 0.5*DEM_MASS;

const double BALL_CENTRE[NDIM] = {0.5,0.5,1.1};
const double BALL_RADIUS = 0.15;

#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY
//#define FIXED_DEM

//#define CHECK_FOR_NAN
//const double POROSITY = 0.90;
const double MAXTIME = 15.0;
const double TIME_DROP_PARTICLE = 5.0;
const int DAMPTIME = 1;
const int OUTSTEP = 500;
const int RESTART_EVERY_N_STEPS = 100;
const int REINIT_DENS_EVERY_N_STEPS = 500000000;



//const double Ar = 8*pow(DEM_RADIUS,3)*(DEM_DENS-DENS)*9.81/VISCOSITY

const double REYNOLDS_NUMBER = 0.001;
const double Ar = REYNOLDS_NUMBER*18;
const double VISCOSITY = sqrt(8*pow(DEM_RADIUS,3)*DENS*(DEM_DENS-DENS)*9.81/Ar)/DENS;
const double VREF = 0.5*VISCOSITY*REYNOLDS_NUMBER/DEM_RADIUS;
const double VMAX = sqrt(2*9.81*(1.0+0.5*POROSITY)*(RMAX[2]-RMIN[2]));

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;
const double SPHBOUNDARYDENS = REFD*pow(9.81*REFD*(RMAX[2]-RMIN[2])/PRB + 1,1.0/7.0);


const double MY_VAR_RES_DV = VMAX/8.0;

const int NX = 30;
const int NY = NX;
const int NZ = NX*1.5;
const double PSEP = (BMAX[0]-BMIN[0])/(NX);
const int MAX_NUM_PARTICLES_PER_CPU = 1.2*NX*NY*NZ;
const double BFAC = 1.00;
//const double H = HFAC*PSEP;
const double H = HFAC*PSEP;


const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


#define MIN_ALPHA 0.1
#endif
