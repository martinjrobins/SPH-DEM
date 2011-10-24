#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

const int QNUMBOTTOM = int(floor(H1/PSEP+0.5));
const int QNUMTOP = int(floor(H2/PSEP+0.5));
const double BSEPBOTTOM = H1/QNUMBOTTOM;
const double BSEPTOP = H2/QNUMTOP;
const double QSEP = 1*PSEP;
const double QSEP2 = 0.5*PSEP;

double qTopLeft[100];
double qBottomLeft[100];
double qBottomRight[100];
double qTopRight[100];



void setBoundary(Cparticle &p,CglobalVars &g,double v1,double v2,double dt) {
   if (p.iam==boundaryobj2) {
      p.r[0] += dt*v2;

      if ((p.norm1[0] < 0 )||(p.norm2[0] < 0)) { //left side of top wedge
         int index = int(floor((p.r[1]-(BMAX[1]-H2))/BSEPTOP+0.5));
         qTopLeft[index] = tan(THETA)*(p.r[1]-(BMAX[1]-H2)) - p.r[0];
      } else if (p.norm1[0] > 0) { //right side of top wedge
         int index = int(floor((p.r[1]-(BMAX[1]-H2))/BSEPTOP+0.5));
         qTopRight[index] = p.r[0] + tan(THETA)*(p.r[1]-(BMAX[1]-H2)) - BMAX[0];
      }
   
   } else if (p.iam==boundaryobj1) {
      p.r[0] += dt*v1;

      if ((p.norm1[0] < 0 )||(p.norm2[0] < 0)) { //left side of bottom wedge
         int index = int(floor((p.r[1]-BMIN[1])/BSEPBOTTOM+0.5));
         qBottomLeft[index] = -tan(THETA)*(p.r[1]-H1-BMIN[1]) - p.r[0];
      } else if (p.norm1[0] > 0) { //right side of bottom wedge
         int index = int(floor((p.r[1]-BMIN[1])/BSEPBOTTOM+0.5));
         qBottomRight[index] = p.r[0] - tan(THETA)*(p.r[1]-H1-BMIN[1]) - BMAX[0];
      }
   }
}


void setBoundaryNormals(Cparticle &p,CglobalVars &g) {
   if (p.iam==boundary) {
      if (p.r[1] > 0.5*(BMAX[1]-BMIN[0])) {
         if (((p.norm1[0] > 0)&&(p.norm1[1] < 0)) || ((p.norm2[0] > 0)&&(p.norm2[1] < 0))) { //top left angle boundary
            int index = int(floor((p.r[1]-(BMAX[1]-H2))/BSEPTOP+0.5));
            if (index != QNUMTOP) {
               double theta;
               if (abs(qTopLeft[index]) < QSEP) {
                  theta = THETA + (0.5*PI-THETA)*((abs(qTopLeft[index])-QSEP)/(QSEP2-QSEP));
                  if (theta >= 0.5*PI) theta = 0.5*PI - 0.000001;
                  p.norm1 = cos(theta),-sin(theta);
                  p.norm2 = 0,0;
               } else {
                  p.norm1 = cos(THETA),-sin(THETA);
                  if (index==0) {
                     p.norm2 = 1,0;
                  } else {
                     p.norm2 = 0,0;
                  }
               }
            }
         } else if (((p.norm1[0] < 0)&&(p.norm1[1] < 0)) || ((p.norm2[0] < 0)&&(p.norm2[1] < 0))) { //top right angle boundary
            int index = int(floor((p.r[1]-(BMAX[1]-H2))/BSEPTOP+0.5));
            if (index != QNUMTOP) {
               double theta;
               if (abs(qTopRight[index]) < QSEP) {
                  theta = THETA + (0.5*PI-THETA)*((abs(qTopRight[index])-QSEP)/(QSEP2-QSEP));
                  if (theta >= 0.5*PI) theta = 0.5*PI - 0.000001;
                  p.norm1 = -cos(theta),-sin(theta);
                  p.norm2 = 0,0;
               } else {
                  p.norm1 = -cos(THETA),-sin(THETA);
                  if (index==0) {
                     p.norm2 = -1,0;
                  } else {
                     p.norm2 = 0,0;
                  }
               }
            }
         }
      } else {
         if (((p.norm1[0] > 0)&&(p.norm1[1] > 0)) || ((p.norm2[0] > 0)&&(p.norm2[1] > 0))) { //bottom left angle boundary
            int index = int(floor((p.r[1]-BMIN[1])/BSEPBOTTOM+0.5));
            if (index != 0) {
               double theta;
               if (abs(qBottomLeft[index]) < QSEP) {
                  theta = THETA + (0.5*PI-THETA)*((abs(qBottomLeft[index])-QSEP)/QSEP2-QSEP);
                  if (theta >= 0.5*PI) theta = 0.5*PI - 0.000001;
                  p.norm1 = cos(theta),sin(theta);
                  p.norm2 = 0,0;
               } else {
                  p.norm1 = cos(THETA),sin(THETA);
                  if (index==QNUMBOTTOM) {
                     p.norm2 = 1,0;
                  } else {
                     p.norm2 = 0,0;
                  }
               }
            }

         } else if (((p.norm1[0] < 0)&&(p.norm1[1] > 0)) || ((p.norm2[0] < 0)&&(p.norm2[1] > 0))) { //bottom right angle boundary
            int index = int(floor((p.r[1]-BMIN[1])/BSEPBOTTOM+0.5));
            if (index != 0) {
               double theta;
               if (abs(qBottomRight[index]) < QSEP) {
                  theta = THETA + (0.5*PI-THETA)*((abs(qBottomRight[index])-QSEP)/(QSEP2-QSEP));
                  if (theta >= 0.5*PI) theta = 0.5*PI - 0.000001;
                  p.norm1 = -cos(theta),sin(theta);
                  p.norm2 = 0,0;
               } else {
                  p.norm1 = -cos(THETA),sin(THETA);
                  if (index==QNUMBOTTOM) {
                     p.norm2 = -1,0;
                  } else {
                     p.norm2 = 0,0;
                  }
               }
            }
         }
      }
   } else if (p.iam==boundaryobj2) {
      if ((p.norm1[0] < 0 )||(p.norm2[0] < 0)) { //left side of top wedge
         int index = int(floor((p.r[1]-(BMAX[1]-H2))/BSEPTOP+0.5));
         if (index != QNUMTOP) {
            double theta;
            if (abs(qTopLeft[index]) < QSEP) {
               theta = THETA + (0.5*PI-THETA)*((abs(qTopLeft[index])-QSEP)/(QSEP2-QSEP));
               if (theta >= 0.5*PI) theta = 0.5*PI - 0.000001;
               p.norm1 = -cos(theta),-sin(theta);
               p.norm2 = 0,0;
            } else {
               p.norm1 = -cos(THETA),-sin(THETA);
               if (index==0) {
                  p.norm2 = 0,-1;
               } else {
                  p.norm2 = 0,0;
               }
            }
         }
      } else if (p.norm1[0] > 0) { //right side of top wedge
         int index = int(floor((p.r[1]-(BMAX[1]-H2))/BSEPTOP+0.5));
         if (index != QNUMTOP) {
            double theta;
            if (abs(qTopRight[index]) < QSEP) {
               theta = THETA + (0.5*PI-THETA)*((abs(qTopRight[index])-QSEP)/(QSEP2-QSEP));
               if (theta >= 0.5*PI) theta = 0.5*PI - 0.000001;
               p.norm1 = cos(theta),-sin(theta);
               p.norm2 = 0,0;
            } else {
               p.norm1 = cos(THETA),-sin(THETA);
               if (index==0) {
                  p.norm2 = 0,-1;
               } else {
                  p.norm2 = 0,0;
               }
            }
         }
      }
   } else if (p.iam==boundaryobj1) {

      if ((p.norm1[0] < 0 )||(p.norm2[0] < 0)) { //left side of bottom wedge
         int index = int(floor((p.r[1]-BMIN[1])/BSEPBOTTOM+0.5));
         if (index != 0) {
            double theta;
            if (abs(qBottomLeft[index]) < QSEP) {
               theta = THETA + (0.5*PI-THETA)*((abs(qBottomLeft[index])-QSEP)/(QSEP2-QSEP));
               if (theta >= 0.5*PI) theta = 0.5*PI - 0.000001;
               p.norm1 = -cos(theta),sin(theta);
               p.norm2 = 0,0;
            } else {
               p.norm1 = -cos(THETA),sin(THETA);
               if (index==QNUMBOTTOM) {
                  p.norm2 = 0,1;
               } else {
                  p.norm2 = 0,0;
               }
            }
         }
      } else if (p.norm1[0] > 0) { //right side of bottom wedge
         int index = int(floor((p.r[1]-BMIN[1])/BSEPBOTTOM+0.5));
         if (index != 0) {
            double theta;
            if (abs(qBottomRight[index]) < QSEP) {
               theta = THETA + (0.5*PI-THETA)*((abs(qBottomRight[index])-QSEP)/(QSEP2-QSEP));
               if (theta >= 0.5*PI) theta = 0.5*PI - 0.000001;
               p.norm1 = cos(theta),sin(theta);
               p.norm2 = 0,0;
            } else {
               p.norm1 = cos(THETA),sin(THETA);
               if (index==QNUMBOTTOM) {
                  p.norm2 = 0,1;
               } else {
                  p.norm2 = 0,0;
               }
            }
         }
      }
   }
}


void CcustomSim::updateSim(const double newTime) {
   double dt = newTime-time;
   double v1,v2;
   if (time<SPEEDUP+DAMP) {
      if (time<DAMP) {
         v1 = 0;
         v2 = 0;
      } else {
         v1 = V1*(time-DAMP)/SPEEDUP;
         v2 = V2*(time-DAMP)/SPEEDUP;
      }
   } else {
      v1 = V1;
      v2 = V2;
   }

   data->traverse<setBoundary,Nsph::ifBoundary>(v1,v2,dt);
   data->traverse<setBoundaryNormals,Nsph::ifBoundary>();
   time = newTime;
}


void CcustomSim::beforeMiddle(double newTime) {
   updateSim(newTime);
}

void CcustomSim::beforeEnd(double newTime) {
   updateSim(newTime);
}

