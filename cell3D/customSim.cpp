#include "customSim.h"
#include "sph.h"
#include "dataLL.impl.h"

inline void addGravity(Cparticle &p,CglobalVars &g) {
   p.ff[2] -= 9.81;
   p.f[2] -= 9.81;
}
inline void setupInflow(Cparticle &p,CglobalVars &g) {
   p.v = 0,0,INFLOW_VEL;
   p.vhat = 0,0,INFLOW_VEL;
}
inline void setH(Cparticle &p,CglobalVars &g) {
   p.h = H;
}
#ifdef LIQ_DEM
inline void addCustomBoundaries(Cparticle &p,CglobalVars &g) {
   const double overlap1 = CYLINDER_ORIGIN[2]+DEM_RADIUS-p.r[2];
   const double overlap2 = p.r[2]-CYLINDER_ORIGIN[2]-CYLINDER_HEIGHT+DEM_RADIUS;
   const double overlap3 = sqrt(pow(p.r[0],2)+pow(p.r[1],2))-CYLINDER_RADIUS+DEM_RADIUS;
   if (overlap1>0) {
      const double overlap_dot = -p.vhat[2];
      vect normal = 0.0;
      normal[2] = 1.0;
      p.fb += Nsph::demNormal(overlap1,overlap_dot,normal)/p.mass;
      //cout << "overlap = "<<overlap1<<" overlap_dot = "<< overlap_dot<<" fb = "<<p.fb<<endl;
   } else if (overlap2>0) {
      const double overlap_dot = p.vhat[2];
      vect normal = 0.0;
      normal[2] = -1.0;
      p.fb += Nsph::demNormal(overlap2,overlap_dot,normal)/p.mass;
   }
   if (overlap3>0) {
      const double theta = atan2(p.r[1],p.r[0]);
      vect normal = 0.0;
      normal[0] = -cos(theta);
      normal[1] = -sin(theta);
      const double overlap_dot = -dot(p.vhat,normal);
      p.fb += Nsph::demNormal(overlap3,overlap_dot,normal)/p.mass;
   }
   p.f += p.fb;

}
#endif

inline void setupInflowTransition(Cparticle &p,CglobalVars &g) {
   const double scaledr = (p.r[2]-(CYLINDER_ORIGIN[2]-INLET_HEIGHT+3.0*PSEP))/(3.0*PSEP);
   const double scaledr2 = (p.r[2]-(CYLINDER_ORIGIN[2]+CYLINDER_HEIGHT+INLET_HEIGHT))/(3.0*PSEP);
   if ((scaledr > 0.0) && (scaledr < 1.0)) {
      p.iam = sph;
      const vect inflowV(0,0,INFLOW_VEL);
      const double ratio = 0.5*(1+cos(PI*scaledr));
      p.v = ratio*inflowV + (1-ratio)*p.v;
      p.vhat = p.v;
   } else if (scaledr2>0.0) {
      p.tag = -111;

      /*
   } else if ((scaledr2 > 0.0) && (scaledr2 < 1.0)) {
      p.iam = sph;
      //get the correct outflow velocity!
      const vect outflowV(0,0,INFLOW_VEL/2.0);   
      const double ratio = 0.5*(1+cos(PI*scaledr2));
      p.v = (1-ratio)*outflowV + ratio*p.v;
      p.vhat = p.v;
   } else if (scaledr2>=1.0) {
      p.iam = sphBoundary;
      const vect outflowV(0,0,INFLOW_VEL/2.0);   
      p.v = outflowV;
      p.vhat = outflowV;
      */
   }
}

void sphBoundaryCircle(CdataLL *data,const vect &origin,const double radiusMin,const double radiusMax) {
   Cparticle p;
   const double dRadius = radiusMax-radiusMin;
   const double NR = floor(dRadius/PSEP+0.5);
   const double RPSEP = dRadius/NR;
   //cout <<"creating inflow with radius = "<<dRadius<<" number of circles = "<<NR<<" separation between circles = "<<RPSEP<<endl;
   const vect n(1,0,0);
   for (int i=0;i<=NR;i++) {
      const double r = i*RPSEP+radiusMin;
      double nc = ceil(2*PI*r/PSEP);
      if (nc == 0) nc++;
      const double TPSEP = 2.0*PI/nc;
      for (int j=0;j<nc;j++) {
         const double theta = j*TPSEP;
         vect newN;
         newN[0] = n[0]*cos(theta)-n[1]*sin(theta);
         newN[1] = n[0]*sin(theta)+n[1]*cos(theta);
         newN[2] = 0;
         p.r = r*newN + origin;
         p.tag = data->globals.time;
         p.mass = pow(PSEP,NDIM)*DENS;
         p.dens = p.mass/pow(RPSEP,NDIM);
         p.h = H;
         p.v = 0,0,0;
         p.iam = sphBoundary;
         data->insertNewParticle(p);
      }
   }
}

inline void driftR(Cparticle &p,CglobalVars &g,double dt) {
   p.r += dt*p.v;
}

void CcustomSim::afterEnd(double newTime) {
   if ((newTime>TIME_START_INLET)&&(floor(time/(PSEP/INFLOW_VEL))<floor(newTime/(PSEP/INFLOW_VEL)))) {
      if (start) {
         vect newOrigin = CYLINDER_ORIGIN;
         newOrigin[2] -= INLET_HEIGHT;
         sphBoundaryCircle(data,newOrigin,0,int((INLET_RADIUS-1.5*PSEP)/PSEP)*PSEP);
      }
      start = true;
      data->traverse<setupInflow,Nsph::ifSphBoundary>();
   }
   time = newTime;
}

void CcustomSim::beforeMiddle(double newTime) {
   if (~doOnce) {
      data->traverse<setH,Nsph::ifDem>();
      doOnce = true;
   }
   data->globals.maxdt = 0.5*PSEP/INFLOW_VEL;
   data->traverse<addGravity,Nsph::ifSphOrDem>();
#ifdef LIQ_DEM
   data->traverse<addCustomBoundaries,Nsph::ifDem>();
#endif
   double dt = newTime-time;
   if ((newTime>TIME_START_INLET)&&(floor(time/(PSEP/INFLOW_VEL))<floor(newTime/(PSEP/INFLOW_VEL)))) {
      if (start) {
         vect newOrigin = CYLINDER_ORIGIN;
         newOrigin[2] -= INLET_HEIGHT;
         sphBoundaryCircle(data,newOrigin,0,int((INLET_RADIUS-1.5*PSEP)/PSEP)*PSEP);
      }
      start = true;
      data->traverse<setupInflow,Nsph::ifSphBoundary>();
   }
   time = newTime;
   //data->traverse<driftR,Nsph::ifSphBoundary>(dt);
}
void CcustomSim::beforeEnd(double newTime) {
   data->traverse<setupInflowTransition,Nsph::ifSphOrSphBoundary>();
   double dt = newTime-time;
   if ((newTime>TIME_START_INLET)&&(floor(time/(PSEP/INFLOW_VEL))<floor(newTime/(PSEP/INFLOW_VEL)))) {
      if (start) {
         vect newOrigin = CYLINDER_ORIGIN;
         newOrigin[2] -= INLET_HEIGHT;
         sphBoundaryCircle(data,newOrigin,0,int((INLET_RADIUS-1.5*PSEP)/PSEP)*PSEP);
      }
      start = true;
      data->traverse<setupInflow,Nsph::ifSphBoundary>();
   }
   time = newTime;
   //data->traverse<driftR,Nsph::ifSphBoundary>(dt);
}
