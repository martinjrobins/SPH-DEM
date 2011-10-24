
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;
const double DENS = 1000.0;
const double DEM_DENS = 2500.0;
const double REFD = 990.0;
#include "parameters.h" 
const double VISCOSITY = 8.92635148e-4/DENS;
const double DEM_RADIUS = 0.0001;
const double AAA = PI*pow(DEM_RADIUS,2.0);
const double DVDV = REYNOLDS_NUMBER*VISCOSITY/(2.0*DEM_RADIUS);
const double DRAGGAMMADRAGGAMMA = 3.7 - 0.65*exp(-pow(1.5-log10(REYNOLDS_NUMBER),2)/2.0);
const double CCC = pow(0.63 + 4.8/pow(REYNOLDS_NUMBER,0.5),2);
const double FDFD = (1.0/2.0)*AAA*CCC*DVDV*DVDV*pow(POROSITY,-DRAGGAMMADRAGGAMMA);

const double RMIN[NDIM] = {0,0,0};
const double RMAX[NDIM] = {0.0066,0.0066,0.0066};
const double BMIN[NDIM] = {0,0,0};
const double BMAX[NDIM] = {0.0066,0.0066,0.0066};


#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 1.4;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {1,1,1};
const int GHOST[2*NDIM] = {0,0,0,0,0,0};

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
//#define MANY_PARTICLES
//#define DENS_DIFFUSE


const double DEM_K = 1.0*10e-5;
const double DEM_VOL = (4.0/3.0)*PI*pow(DEM_RADIUS,3);
const double DEM_MASS = DEM_VOL*DEM_DENS;
const double DEM_GAMMA = 2.0*10e-9;
const double DEM_MIN_REDUCED_MASS = 0.5*DEM_MASS;
#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY
#define FIXED_DEM

//#define CHECK_FOR_NAN
const double MAXTIME = (RMAX[0]-RMIN[0])/(2.0*DVDV);
const double TIME_DROP_PARTICLE = 0.0;
const int DAMPTIME = 0.0;
const int OUTSTEP = 200;
const int RESTART_EVERY_N_STEPS = 250;
const int REINIT_DENS_EVERY_N_STEPS = 500000000;


//const double Ar = 8*pow(DEM_RADIUS,3)*(DEM_DENS-DENS)*9.81/VISCOSITY

const double VREF = 10.0*DVDV;
const double VMAX = 10.0*DVDV;

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;
const double SPHBOUNDARYDENS = REFD*pow(9.81*REFD*(RMAX[2]-RMIN[2])/PRB + 1,1.0/7.0);

const double DPRESS = ((1-POROSITY)/DEM_VOL)*FDFD*(RMAX[0]-RMIN[0]);
const double D2 = REFD*pow((DPRESS/PRB)+pow(DENS/REFD,GAMMA),1.0/GAMMA);
const double DENS_DROP[NDIM] = {POROSITY*(DENS-D2),0,0};


const double MY_VAR_RES_DV = VMAX/8.0;

const int NY = NX;
const int NZ = NX;
const double PSEP = (BMAX[0]-BMIN[0])/(NX);
const int MAX_NUM_PARTICLES_PER_CPU = NX*NY*NZ+1.1*(RMAX[0]-RMIN[0])*(RMAX[1]-RMIN[1])*(RMAX[2]-RMIN[2])*(1-POROSITY)/DEM_VOL;
const double BFAC = 1.00;
//const double H = HFAC*PSEP;
const double H = HFAC*PSEP;


const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


#define MIN_ALPHA 0.1
#endif
