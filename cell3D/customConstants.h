
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>
#include <algorithm>

#define _3D_
#define NDIM 3
#define INCOMPRESS

const double PI = 3.14159265;

#include "parameters.h"

#if NIX == 10
const double INLET_N = 98;
#elif NIX == 8
const double INLET_N = 66;
#elif NIX == 6
const double INLET_N = 40;
#elif NIX == 4
const double INLET_N = 21;
#elif NIX == 3
const double INLET_N = 14;
#elif NIX == 2
const double INLET_N = 8;
#else
const double INLET_N = 1;
#endif

const double CYLINDER_ORIGIN[NDIM] = {0.0,0.0,0.0};
const double CYLINDER_RADIUS = 11.5/1000.0;
const double CYLINDER_HEIGHT = 31.0/1000.0;
const double TOP_BUFFER = 0.05;
//const double INLET_RADIUS = 1.9/1000.0;
const double BOUNDARY_BUFFER = 1.1;
const double PSEP = REAL_INLET_RADIUS*2.0/(NIX+BOUNDARY_BUFFER);
const double INLET_RADIUS = REAL_INLET_RADIUS+0.5*BOUNDARY_BUFFER*PSEP;
//const double REAL_INLET_RADIUS = 0.5/1000.0;
const double INLET_HEIGHT = 6*PSEP;
//const double INLET_FLOW_RATE = 100.0;
const double INFLOW_VEL = INLET_FLOW_RATE*1e-6/(PI*pow(REAL_INLET_RADIUS,2)*60.0);
const double INLET_AREA = INLET_N*PSEP*PSEP;
//const double INFLOW_VEL = 1e-6*(1.0/60.0)*INLET_FLOW_RATE/INLET_AREA;
const double RMIN[NDIM] = {CYLINDER_ORIGIN[0]-CYLINDER_RADIUS,CYLINDER_ORIGIN[1]-CYLINDER_RADIUS,CYLINDER_ORIGIN[2]-INLET_HEIGHT};
const double RMAX[NDIM] = {CYLINDER_ORIGIN[0]+CYLINDER_RADIUS,CYLINDER_ORIGIN[1]+CYLINDER_RADIUS,CYLINDER_ORIGIN[2]+CYLINDER_HEIGHT*(1+TOP_BUFFER)+INLET_HEIGHT};
const double BMIN[NDIM] = {CYLINDER_ORIGIN[0]-CYLINDER_RADIUS,CYLINDER_ORIGIN[1]-CYLINDER_RADIUS,CYLINDER_ORIGIN[2]-INLET_HEIGHT};
const double BMAX[NDIM] = {CYLINDER_ORIGIN[0]+CYLINDER_RADIUS,CYLINDER_ORIGIN[1]+CYLINDER_RADIUS,CYLINDER_ORIGIN[2]+CYLINDER_HEIGHT*(1+TOP_BUFFER)+INLET_HEIGHT};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

const double GAMMA = 7.0;
const double HFAC = 1.5;

const int PERIODIC[NDIM] = {0,0,0};
const int GHOST[2*NDIM] = {0,0,0,0,0,0};
const double DENS_DROP[NDIM] = {0,0,0};

const double DENS = 997.0479;
const double REFD = 997.0479;
//#define QUINTIC
#define WENDLAND
#define VISC_MORRIS
#define VISC_ARTIFICIAL
const double ALPHA_ARTIFICIAL = 0.1;
//#define HANN
//#define VISC_MONAGHAN
//#define VISC_CLEARY
//#define VORT_LEASTSQUARES
#define VAR_H_CORRECTION
//#define REFERENCE_FRAME
//

#define LIQ_DEM
#define LIQ_DEM_CUSTOM_WALL_CONTACTS
#define LIQ_DEM_SEPARATE_DRAGS
#define LIQ_DEM_VARIABLE_TIMESTEP
//#define LINEAR
//#define LUBRICATION
const double ROUGHNESS_EPSILON = 0.01;
#define NO_ANTICLUMPING
//#define CONST_H
//#define MORRIS_SPH_BOUNDARY

//#define BOUNDARY_FORCE_WITH_HEIGHT

//#define CHECK_FOR_NAN
const double MAXTIME = 8.0;
const double DAMPTIME = 0.1;
const int OUTSTEP = 1000;
const int RESTART_EVERY_N_STEPS = 99;
const int REINIT_DENS_EVERY_N_STEPS = 2000000;
const int REINIT_DENS_AT_N_STEPS = 2000000;
const double TIME_START_INLET = MAXTIME/10.0;
const double TIME_CHANGE_DEM_MASS = 4.0*TIME_START_INLET/5.0;

const double VREF = INFLOW_VEL;
//const double VREF = 1.0;
const double VMAX = std::max(VREF,sqrt(2.0*9.81*(CYLINDER_HEIGHT+INLET_HEIGHT)));
const double VISCOSITY = 8.92635148e-4/DENS;
const double REYNOLDS_NUMBER = VREF*(BMAX[0]-BMIN[0])/VISCOSITY;

//const double LYAP_DT = 1.0;

const double CSFAC = 10.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,GAMMA-1.0)*pow(SPSOUND,2)*REFD/GAMMA;


const double MY_VAR_RES_DV = VMAX/8.0;

//const int NX = 30;
const int NX = (BMAX[0]-BMIN[0])/PSEP;
const int NY = NX;
const int NZ = CYLINDER_HEIGHT/PSEP;
const double BFAC = 1.00;
//const double H = HFAC*PSEP;
const double H = HFAC*PSEP;

//const double DEM_RADIUS = (BMAX[0]-BMIN[0])/200;
//const double DEM_RADIUS = 0.167/2.0;
//const double DEM_RADIUS = 0.65/1000.0;
//const double DEM_DENS = DENS*1.5;
const double DEM_MASS = (4.0/3.0)*PI*pow(DEM_RADIUS,3)*DEM_DENS;
const double INIT_DEM_DENS = DEM_DENS*3.0;
const double INIT_DEM_MASS = DEM_MASS*3.0;
const double DEM_DISK_HEIGHT = CYLINDER_HEIGHT/20.0;
const double DEM_DISK_POROSITY = 0.5;
//const double DEM_K = 1.0e03;
//const double DEM_GAMMA = 0.004;
const double DEM_K = 1.0e01;
const double DEM_GAMMA = 0.0004;

const double DEM_MIN_REDUCED_MASS = 0.5*DEM_MASS;
//const double DEM_COEFF_REST = exp(-PI*DEM_GAMMA/(2.0*DEM_MIN_REDUCED_MASS*sqrt(DEM_K/DEM_MIN_REDUCED_MASS-pow(DEM_GAMMA/(2.0*DEM_MIN_REDUCED_MASS),2))));
//const double DEM_GAMMA = 2.0*DEM_MIN_REDUCED_MASS*sqrt(DEM_K/DEM_MIN_REDUCED_MASS - pow(10000*PI/20.0,2));
//const double DEM_BED_HEIGHT = 3.0/1000;
//const double DEM_BED_POROSITY = 0.37;
//const double SOLIDVOL = PI*pow(CYLINDER_RADIUS,2)*DEM_BED_HEIGHT*(1-DEM_BED_POROSITY);
const double DEM_TOTAL_MASS = 6.0/1000.0;
const double DEM_VOL = (4.0/3.0)*PI*pow(DEM_RADIUS,3);
#ifdef LIQ_DEM
const int NDEM = DEM_TOTAL_MASS/(DEM_VOL*DEM_DENS);
#else
const int NDEM = 0;
#endif
const double DEM_DISK_VOL = PI*pow(CYLINDER_RADIUS,2)*DEM_DISK_HEIGHT;
const double DEM_DISK_NDEM = DEM_DISK_VOL*(1-DEM_DISK_POROSITY)/DEM_DISK_VOL;
const double TIME_DROP_PARTICLE = 0.0;

const double Konmr = DEM_K/DEM_MIN_REDUCED_MASS;
const double gammaOnmr = pow(0.5*DEM_GAMMA/DEM_MIN_REDUCED_MASS,2);
const double THETS = PI/sqrt(Konmr-gammaOnmr);
const double CR = exp(-DEM_GAMMA*THETS/(2.0*DEM_MIN_REDUCED_MASS));

const int NCPU = NCPU_X*NCPU_Y*NCPU_Z;
const int MAX_NUM_PARTICLES_PER_CPU = 1.3*PI*pow(CYLINDER_RADIUS,2)*(CYLINDER_HEIGHT*(1+TOP_BUFFER))/(pow(PSEP,NDIM)*NCPU) + 1.0*NDEM;
const int MAX_NUM_DEM_PARTICLES = 1.0*NDEM;
//const double LIQ_DEM_COUPLING_RADIUS = 6.8*DEM_RADIUS;
const double LIQ_DEM_COUPLING_RADIUS = 8*DEM_RADIUS;

const double APPROX_SMOOTH_GRID_SIZE = H;
const double SMOOTH_GRID_SIZE = (BMAX[0]-BMIN[0])/floor((BMAX[0]-BMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;


#define MIN_ALPHA 0.1
#endif
