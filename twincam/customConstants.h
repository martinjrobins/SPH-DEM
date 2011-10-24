
#ifndef CUSTOMCONSTANTS_H
#define CUSTOMCONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS


const double RMIN[NDIM] = {0.0,0.0};
const double RMAX[NDIM] = {1.2,0.6};

#define TIME_OVERIDE 0
#define START_TIME 0
#define NSTEP_OVERIDE 0
#define START_NSTEP 0

#define GAMMA 1.4
#define HFAC  1.95
#define PI    3.14159265

const int PERIODIC[NDIM] = {0,0};

#define DENS 1500.0
#define REFD 1500.0
//#define QUINTIC
//#define VISC_MORRIS
#define VISC_MONAGHAN
//#define VISC_CLEARY
#define VORT_LEASTSQUARES
//#define REFERENCE_FRAME

const double BOX_GLOB_TIME = 0;
const double PERIOD = 120;
const double NOREVS = 3;
const double MAXTIME = BOX_GLOB_TIME+ NOREVS*PERIOD;
const int OUTSTEP = 400;
const int RESTART_EVERY_N_STEPS = 100;
const int REINIT_DENS_EVERY_N_STEPS = 200000000;

const double S = 0.45;
const double TANK_H_OVER_WIDTH = 60.0/108.0;
const double TANK_WIDTH = RMAX[0]-RMIN[0];
const double TANK_R = TANK_WIDTH*TANK_H_OVER_WIDTH/2.0;
const double TANK_C1[2] = {TANK_R+RMIN[0],TANK_R+RMIN[1]};
const double TANK_C2[2] = {TANK_WIDTH+RMIN[0]-TANK_R,TANK_R+RMIN[1]};

const double CAM_R = 1.0/sqrt(3.0)*S;
const double CAM1_RPM = 0.5;
const double CAM2_RPM = 0.5;
const double CAM1_ANGVEL = CAM1_RPM/60.0*2.0*PI;
const double CAM2_ANGVEL = CAM1_RPM/60.0*2.0*PI;
const double CAM1_PHASE = 0;
const double CAM2_PHASE = 0;
const double CAM1_ANG = 60;
const double CAM2_ANG = 60;
const double VREF = CAM1_ANGVEL*CAM_R;
const double VMAX = VREF;

const double VISCOSITY = 0.005;
const double REYNOLDS_NUMBER = VREF*CAM_R/VISCOSITY;

const double CSFAC = 10.0;
const double SPSOUND = CSFAC*VMAX;
const double PRB = pow(REFD/DENS,6)*pow(SPSOUND,2)*REFD/7.0;
//const double PRB = pow(SPSOUND,2)*REFD/7.0;


//const double ALPHA = 0.1;
//const double H = Nmisc::viscToH(VISCOSITY,ALPHA,SPSOUND);
//const double PSEP = H/HFAC;
//const int NX = static_cast<int>((RMAX[0]-RMIN[0])/PSEP);
//const int NY = static_cast<int>((RMAX[1]-RMIN[1])/PSEP);

const int NX = 160;
const double PSEP = (RMAX[0]-RMIN[0])/NX;
//const double H = HFAC*PSEP;
const double H = HFAC*(RMAX[0]-RMIN[0])/NX;

const int NY = int(TANK_R*2.0/PSEP);

const double APPROX_SMOOTH_GRID_SIZE = H/2;
const double SMOOTH_GRID_SIZE = (RMAX[0]-RMIN[0])/floor((RMAX[0]-RMIN[0])/APPROX_SMOOTH_GRID_SIZE+0.5);

const double GRIDSEP = PSEP;

#define MIN_ALPHA 0.1
#endif
