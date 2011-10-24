
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;

const double CYLINDER_ORIGIN[NDIM] = {0.0,0.0,0.0};
const double CYLINDER_RADIUS = 11.5;
const double CYLINDER_HEIGHT = 31.0;
const double RMIN[NDIM] = {-CYLINDER_RADIUS-0.1,-CYLINDER_RADIUS-0.1,-0.1};
const double RMAX[NDIM] = {CYLINDER_RADIUS+0.1,CYLINDER_RADIUS+0.1,CYLINDER_HEIGHT+0.1};
const double BMIN[NDIM] = {-CYLINDER_RADIUS,-CYLINDER_RADIUS,0.0};
const double BMAX[NDIM] = {CYLINDER_RADIUS,CYLINDER_RADIUS,CYLINDER_HEIGHT};

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
//

#define LIQ_DEM
#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY


const double MAXTIME = 20;
const int DAMP_STEPS = 100;
const int OUTSTEP = 500;
const int RESTART_EVERY_N_STEPS = 100;
const int REINIT_DENS_EVERY_N_STEPS = 2000000;

const double VREF = 10.0;
const double VMAX = VREF;
const double VISCOSITY = 0.892635148;
const double REYNOLDS_NUMBER = VREF*(BMAX[0]-BMIN[0])/VISCOSITY;

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;


const double MY_VAR_RES_DV = VMAX/8.0;

const int NX = 10;
const int NY = NX;
const double PSEP = 0.5*(BMAX[0]-BMIN[0])/(NX+0.3);
const int NZ = CYLINDER_HEIGHT*0.7/PSEP;
const double BFAC = 1.00;
//const double H = HFAC*PSEP;
const double H = HFAC*PSEP;

//const double DEM_RADIUS = (BMAX[0]-BMIN[0])/200;
const double DEM_RADIUS = 0.167/2.0;
const double DEM_DENS = DENS*10.0;

const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
