#ifndef DEMCONSTANTS_H
#define DEMCONSTANTS_H

#ifdef LIQ_DEM
#ifdef LINEAR
const double DEM_TIMESTEP = (1.0/50.0)*PI*sqrt(DEM_MIN_REDUCED_MASS)/sqrt(DEM_K-pow(0.5*DEM_GAMMA,2)/DEM_MIN_REDUCED_MASS);
#endif
const double LIQ_DEM_TIMESTEP = (1.0/20.0)*DEM_VOL*DEM_DENS/(6.0*PI*VISCOSITY*DENS*DEM_RADIUS);
#endif


//const int TWIN_N = LIQ_DEM_TWIN/DEM_TIMESTEP;
const int TWIN_N = 2;
 
#endif
