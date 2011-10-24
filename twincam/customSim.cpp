#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

bool ifCamBoundary(Cparticle &p) {
   return (p.iam == boundary)&&((pow(p.r[0]-TANK_C1[0],2)+pow(p.r[1]-TANK_C1[1],2) < TANK_R*TANK_R)||(pow(p.r[0]-TANK_C2[0],2)+pow(p.r[1]-TANK_C2[1],2) < TANK_R*TANK_R));
}
bool ifLeftCamBoundary(Cparticle &p) {
   return (p.iam == boundary)&&(p.r[0] < 0.5*(RMIN[0]+RMAX[0]))&&(pow(p.r[0]-TANK_C1[0],2)+pow(p.r[1]-TANK_C1[1],2) < TANK_R*TANK_R);
}
bool ifRightCamBoundary(Cparticle &p) {
   return (p.iam == boundary)&&(p.r[0] > 0.5*(RMIN[0]+RMAX[0]))&&(pow(p.r[0]-TANK_C2[0],2)+pow(p.r[1]-TANK_C2[1],2) < TANK_R*TANK_R);
}

void setBoundary(Cparticle &p,CglobalVars &g,double dt) {
   double rotang,angvel;

   if (ifLeftCamBoundary(p)) {
      angvel = CAM1_ANGVEL;
   } else {
      angvel = CAM2_ANGVEL;
   }
   rotang = angvel*dt;
      
   double cang = cos(rotang);
   double sang = sin(rotang);
   vect r = p.r;
   p.r[0] = r[0]*cang-r[1]*sang;
   p.r[1] = r[0]*sang+r[1]*cang;

   p.v[0] = -angvel*p.r[1];
   p.v[1] = angvel*p.r[0];
}



void CcustomSim::updateSim(const double newTime) {
   double dt = newTime-time;
   double middle_time = time+dt/2;
   
   if (middle_time >= BOX_GLOB_TIME) {
      data->traverse<setBoundary,ifCamBoundary>(dt);
   }
   time = newTime;
}


void CcustomSim::beforeMiddle(double newTime) {
   updateSim(newTime);
}

void CcustomSim::beforeEnd(double newTime) {
   updateSim(newTime);
}

vect vFilter(vect vin) {
   
}
      //vect rFilter(vect) {}

