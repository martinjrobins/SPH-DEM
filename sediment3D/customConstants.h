
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;
#include "parameters.h" 
const double DEM_DENS = 2500;
const double REFD = DENS;
//#define VISCOSITY (0.05 > REYNOLDS_NUMBER ? 8.92635148e-3/DENS : 8.92635148e-4/DENS)
//const double BETA = 3.7 - 0.65*exp(-(pow(1.5-log10(REYNOLDS_NUMBER),2))/2.0);
//const double Ar =  (3.0/4.0)*pow(POROSITY,-BETA)*(0.3969*pow(REYNOLDS_NUMBER,2) + 6.048*pow(REYNOLDS_NUMBER,1.5) + 23.04*REYNOLDS_NUMBER);
//const double DEM_RADIUS = pow(Ar*pow(VISCOSITY,2)/((DEM_DENS/DENS-1)*9.81),1.0/3.0)/2.0;

const double RMIN[NDIM] = {0,0,0};
const double L = 0.004;
const double RMAX[NDIM] = {L,L,L*HMULT};
const double WMAX = L*HMULT-0.5*L;
const double BMIN[NDIM] = {0,0,0};
const double BMAX[NDIM] = {L,L,L*HMULT};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 7.0;
const double HFAC = 1.5;

//#define GHOST_BOUNDARY
const int PERIODIC[NDIM] = {1,1,0};
#ifdef GHOST_BOUNDARY
const int GHOST[2*NDIM] = {0,0,0,0,2,0};
#else
const int GHOST[2*NDIM] = {0,0,0,0,0,0};
#endif
const double DENS_DROP[NDIM] = {0,0,0};

//#define QUINTIC
#define WENDLAND
#define VISC_MORRIS
#define VISC_ARTIFICIAL
const double ALPHA_ARTIFICIAL = 0.1;
//#define HANN
//#define VISC_MONAGHAN
//#define VISC_CLEARY
//#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME
//#define CORRECTED_GRADIENT
#define VAR_H_CORRECTION

//#define SETTLE
//#define SPHBOUNDARY
#define LIQ_DEM
//#define LIQ_DEM_TEST
//#define LIQ_DEM_SEPARATE_DRAGS
//#define HALF_COURANT
//#define LIQ_DEM_SIMPLE_DRAG
//#define LIQ_DEM_ADDED_MASS
//#define LIQ_DEM_ONE_WAY_COUPLE
#define LIQ_DEM_CUSTOM_WALL_CONTACTS
//#define LIQ_DEM_DDDT_VER2
#define LINEAR
//#define MANY_PARTICLES
//#define DENS_DIFFUSE

const double DEM_RADIUS = 50e-6;
const double DEM_K = 1.0*10e-4;
const double DEM_VOL = (4.0/3.0)*PI*pow(DEM_RADIUS,3);
const double DEM_MASS = DEM_VOL*DEM_DENS;
//const double DEM_GAMMA = 1.0*10e-10;
const double DEM_GAMMA = 0;
const double DEM_MIN_REDUCED_MASS = 0.5*DEM_MASS;
#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY
//#define FIXED_DEM

//#define CHECK_FOR_NAN
//const double POROSITY = 0.90;

//const double Ar = 8*pow(DEM_RADIUS,3)*(DEM_DENS-DENS)*9.81/VISCOSITY

const double VREF = 0.5*VISCOSITY*REYNOLDS_NUMBER/(POROSITY*DEM_RADIUS);
//const double VMAX = sqrt(2*9.81*(1.0+0.5*(1-POROSITY)*(DEM_DENS-DENS)/DENS)*(RMAX[2]-RMIN[2]));
const double VMAX = 2.0*sqrt(2*9.81*(RMAX[2]-RMIN[2]));

#ifdef SHORT_TIME
const double MAXTIME = 0.15*L*HMULT/VREF;
#else
const double MAXTIME = 1.5*L*HMULT/VREF+0.1;
#endif
const double TIME_DROP_PARTICLE = (1.0/3.0)*MAXTIME;
const int DAMPTIME = TIME_DROP_PARTICLE/2.0;
#ifdef MANY_PARTICLES
const int OUTSTEP = 200;
#else
const int OUTSTEP = 1000;
#endif
const int RESTART_EVERY_N_STEPS = 250;
const int REINIT_DENS_EVERY_N_STEPS = 50000000;
const int REINIT_DENS_AT_N_STEPS = 50000000;
//#define REINIT_DENS_MLS



//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
//const double PRB = pow(SPSOUND,2)*DENS/(pow(DENS/REFD,GAMMA)*(GAMMA+1)-1);
const double PRB = pow(REFD/DENS,GAMMA-1.0)*pow(SPSOUND,2)*REFD/GAMMA;
const double SPHBOUNDARYDENS = REFD*pow(9.81*REFD*(RMAX[2]-RMIN[2])/PRB + 1,1.0/7.0);


const double MY_VAR_RES_DV = VMAX/8.0;

const int NY = NX;
const double PSEP = (BMAX[0]-BMIN[0])/(NX);
const int NZ = WMAX/PSEP;
const int MAX_NUM_PARTICLES_PER_CPU = 1.1*(NX*NY*NZ+(1-POROSITY)*0.5*(RMAX[0]-RMIN[0])*(RMAX[1]-RMIN[1])*(RMAX[2]-RMIN[2])/DEM_VOL);
const int MAX_NUM_DEM_PARTICLES = 1.1*(1-POROSITY)*0.5*(RMAX[0]-RMIN[0])*(RMAX[1]-RMIN[1])*(RMAX[2]-RMIN[2])/DEM_VOL+1;
const double BFAC = 1.00;
//const double H = HFAC*PSEP;
const double H = HFAC*PSEP;
const double LIQ_DEM_COUPLING_RADIUS = 1.999*H/pow(POROSITY,1.0/NDIM);


const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


#define MIN_ALPHA 0.1
#endif
