
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

#define _2D_
#define NDIM 2
#define INCOMPRESS
#define MAXTIME 10.0
#define MAXSTEP 700
#define IOSTEP  2

#define NX 50
#define NY 20

const double RMIN[NDIM] = {-0.5,-0.5};
const double RMAX[NDIM] = {0.5,0.5};

#define TIME_OVERIDE 1
#define START_TIME 0
#define NSTEP_OVERIDE 1
#define START_NSTEP 0

#define GAMMA 1.4
#define HFAC  1.3
#define PI    3.14159265

const double PMIN[NDIM] = {-0.5,-100};
const double PMAX[NDIM] = {0.5,100};

#define DENS 1500
#define REFD 1500
#define VREF 1.0
#define ALPHA 0.1
#define CSFAC 10
const double PSEP = (RMAX[0]-RMIN[0])/NX;
const double H = HFAC*PSEP;
const double PRB = pow(CSFAC*VREF,2)*REFD/7.0;
const double SPSOUND = sqrt(7.0*PRB/REFD);
const double WALLSPEED = (RMAX[0]-RMIN[0])/1.0;

#define MIN_ALPHA 0.1
#endif
