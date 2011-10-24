#include "customConstants.h"
#include "particle.h"
#include "io_data_vtk.h"
#include "misc.h"
#include "sph.h"
#include "sphIncompress.h"
#include "customOutput.h"
#include "dataLL.h"
#include <vector>
#include <ctime>
#include <iostream>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


bool otherSide(double m,double b,double x1,double y1,double x2,double y2) {
   double mtmp,btmp,xint;
   bool other = false;
   mtmp = (y2-y1)/(x2-x1);
   btmp = y1 - mtmp*x1;
   if (m!=mtmp) {
      xint = (btmp-b)/(m-mtmp);
      if (((xint>x1)&&(xint<x2))||((xint<x1)&&(xint>x2))) {
         other = true;
      }
   }
   //cout << "m= "<<m<<" b = "<<b<<" x1 = "<<x1<<" y1 = "<<y1<<" x2 = "<<x2<<" y2 = "<<y2<<" other = "<<other<<endl;
   return other;
}


bool inTriangle(double xx,double yy,double cenX,double cenY,double ang,double r,double phase) {
   bool isin;

   double baseAng,angRad,ax,ay,bx,by,cx,cy;
   angRad = ang*PI/180.0;
   baseAng = phase*PI/180.0 + PI/2.0;
 
   ax = r*cos(baseAng) + cenX;
   ay = r*sin(baseAng) + cenY;

   bx = r*cos(baseAng+PI-angRad) + cenX;
   by = r*sin(baseAng+PI-angRad) + cenY;

   cx = r*cos(baseAng+PI+angRad) + cenX;
   cy = r*sin(baseAng+PI+angRad) + cenY;

   double m1,b1,m2,b2,m3,b3;

   isin = false;
   if (bx != ax) {
      m1 = (by-ay)/(bx-ax);
      b1 = ay - m1*ax;
      isin = otherSide(m1,b1,xx,yy,cx,cy);
   } else if (((bx>xx)&&(bx<cx))||((bx<xx)&&(bx>cx))) {
      isin = true;
   }

   if (bx != cx) {
      m2 = (cy-by)/(cx-bx);
      b2 = by - m2*bx;
      isin = isin || otherSide(m2,b2,xx,yy,ax,ay);
   } else if (((bx>xx)&&(bx<ax))||((bx<xx)&&(bx>ax))) {
      isin = true;
   }

   if (cx != ax) {
      m3 = (ay-cy)/(ax-cx);
      b3 = cy - m3*cx;
      isin = isin || otherSide(m3,b3,xx,yy,bx,by);
   } else if (((cx>xx)&&(cx<bx))||((cx<xx)&&(cx>bx))) {
      isin = true;
   }

   return !isin;
}


void genTrianglePoints(vector<Cparticle> &ps,double cenX,double cenY,double ang,double r,double phase) {
   double ax,ay,bx,by,cx,cy;
   double baseAng,angRad;
   angRad = ang*PI/180.0;
   baseAng = phase*PI/180.0 + PI/2.0;
 
   ax = r*cos(baseAng) + cenX;
   ay = r*sin(baseAng) + cenY;

   bx = r*cos(baseAng+PI-angRad) + cenX;
   by = r*sin(baseAng+PI-angRad) + cenY;

   cx = r*cos(baseAng+PI+angRad) + cenX;
   cy = r*sin(baseAng+PI+angRad) + cenY;

   cout << "triangle points at (" <<ax<<","<<ay<<"), ("<<bx<<","<<by<<"), ("<<cx<<","<<cy<<")"<<endl;

   double b1 = sqrt(3.0)/6.0*S;
   double b2 = 1.0/sqrt(3.0)*S;
   double xb[4],yb[4],anxb[4],anyb[4];
   double tmp;

   xb[0] = ax;
   yb[0] = ay;
   anxb[0] = by-ay;
   anyb[0] = ax-bx;
   tmp = sqrt(pow(anxb[0],2)+pow(anyb[0],2));
   anxb[0] = anxb[0]/tmp;
   anyb[0] = anyb[0]/tmp;

   xb[1] = bx;
   yb[1] = by;
   anxb[1] = cy-by;
   anyb[1] = bx-cx;
   tmp = sqrt(pow(anxb[1],2)+pow(anyb[1],2));
   anxb[1] = anxb[1]/tmp;
   anyb[1] = anyb[1]/tmp;

   xb[2] = cx;
   yb[2] = cy;
   anxb[2] = ay-cy;
   anyb[2] = cx-ax;
   tmp = sqrt(pow(anxb[2],2)+pow(anyb[2],2));
   anxb[2] = anxb[2]/tmp;
   anyb[2] = anyb[2]/tmp;

   xb[3] = ax;
   yb[3] = ay;
   anxb[3] = by-ay;
   anyb[3] = ax-bx;
   tmp = sqrt(pow(anxb[3],2)+pow(anyb[3],2));
   anxb[3] = anxb[3]/tmp;
   anyb[3] = anyb[3]/tmp;

   for (int i=0;i<3;i++) {
      double slen = sqrt(pow(xb[i+1]-xb[i],2)+pow(yb[i+1]-yb[i],2));
      double drhatx = (xb[i+1]-xb[i])/slen;
      double drhaty = (yb[i+1]-yb[i])/slen;
      double nsnew = int(slen/PSEP+0.5);
      double dpnew = slen/nsnew;
      for (int j=0;j<nsnew;j++) {
         double xx = xb[i]+drhatx*dpnew*j;
         double yy = yb[i]+drhaty*dpnew*j;

         Cparticle p;
         
         p.tag = ps.size()+1;
         p.r = xx,yy;
         p.dens = DENS;
         p.mass = PSEP*PSEP*DENS;
         p.h = H;
         p.v = 0.0,0.0;
         p.iam = boundary;
         p.norm1 = anxb[i],anyb[i];
         if (j==0) {
            p.concave = 0;
            if (i==0) {
               p.norm2 = anxb[2],anyb[2];
            } else {
               p.norm2 = anxb[i-1],anyb[i-1];
            }
         }
         p.alpha = ALPHA;
         ps.push_back(p);
      }
   }
}

void genTriangles(vector <Cparticle> &ps, int &nLeftCam, int &nRightCam) {
   genTrianglePoints(ps,TANK_C1[0],TANK_C1[1],CAM1_ANG,CAM_R,CAM1_PHASE);
   nLeftCam = ps.size();
   genTrianglePoints(ps,TANK_C2[0],TANK_C2[1],CAM2_ANG,CAM_R,CAM2_PHASE);
   nRightCam = ps.size();
}



int main(int argc, char *argv[]) {
   if (argc != 2) {
      cout << "Usage: setup outfilename" << endl;
      return(-1);
   }
   string filename = argv[1];

   vector<Cparticle> ps;
   Cparticle p;
   cout << "Creating simulation with sides: rmax = ["<<RMAX[0]<<" "<<RMAX[1]<<"] rmin = ["<<RMIN[0]<<" "<<RMIN[1]<<"]"<<endl;
   cout << "Reynolds Number = "<<REYNOLDS_NUMBER<<endl;
   cout << "Density = "<<DENS<<endl;
   cout << "number of particles on side = "<<NX<<endl;
   cout << "alpha = "<<ALPHA<<endl;
   cout << "viscosity = "<<VISCOSITY<<endl;
   cout << "maxtime = "<<MAXTIME<<endl;
   cout <<"damp time = "<<BOX_GLOB_TIME<<endl;
 
   double tsc = Nsph::courantCondition(H,2*SPSOUND);
   double tsv = Nsph::viscDiffusionCondition(H,VISCOSITY);
   cout <<"simulation will take "<<int((MAXTIME/tsc)+1)<<" steps according to Courant condition, "<<int((MAXTIME/tsv)+1)<<" steps according to visc diffusion condition"<<endl;
   
   int nLeftCam,nRightCam;
   genTriangles(ps,nLeftCam,nRightCam);
   
   for (int i=0;i<=NX;i++) {
      for (int j=0;j<=NY;j++) {
         double xx = i*PSEP+RMIN[0];
         double yy = j*PSEP+RMIN[1];

         double r1 = pow((xx-TANK_C1[0]),2) + pow((yy-TANK_C1[1]),2);
         double r2 = pow((xx-TANK_C2[0]),2) + pow((yy-TANK_C2[1]),2);

         if ((r1<TANK_R*TANK_R)||(r2<TANK_R*TANK_R)) {
            if ((!inTriangle(xx,yy,TANK_C1[0],TANK_C1[1],CAM1_ANG,CAM_R+3*PSEP,CAM1_PHASE))&&(!inTriangle(xx,yy,TANK_C2[0],TANK_C2[1],CAM2_ANG,CAM_R+3*PSEP,CAM2_PHASE))) {
               p.tag = ps.size()+1;
               p.r = xx,yy;
               p.dens = DENS;
               p.mass = PSEP*PSEP*DENS;
               p.h = H;
               p.v = 0.0,0.0;
               p.iam = sph;

               p.alpha = ALPHA;
               ps.push_back(p);
            }
         }
      }
   }

   const double STHETA = acos((TANK_C2[0]-TANK_C1[0])/(2.0*TANK_R));
   const double ETHETA = 2.0*PI-STHETA;
   const double NTANKPERIM = int((ETHETA-STHETA)*TANK_R/PSEP);
   const double DTHETA = (ETHETA-STHETA)/NTANKPERIM;
   
   for (int i=1;i<=NTANKPERIM;i++) {
      double xx = TANK_R*cos((i-1)*DTHETA+STHETA)+TANK_C1[0];
      double yy = TANK_R*sin((i-1)*DTHETA+STHETA)+TANK_C1[1];

      p.tag = ps.size()+1;
      p.r = xx,yy;
      p.dens = DENS;
      p.mass = PSEP*PSEP*DENS;
      p.h = H;
      p.v = 0.0,0.0;
      p.iam = boundary;
      p.norm1[0] = -cos((i-1)*DTHETA+STHETA);
      p.norm1[1] = -sin((i-1)*DTHETA+STHETA);
      if (i==NTANKPERIM) {
         p.concave = 0;
         p.norm2[0] = -cos(STHETA+PI);
         p.norm2[1] = -sin(STHETA+PI);
      }
      p.alpha = ALPHA;
      ps.push_back(p);

   }


   for (int i=1;i<=NTANKPERIM;i++) {
      double xx = TANK_R*cos((i-1)*DTHETA+STHETA+PI)+TANK_C2[0];
      double yy = TANK_R*sin((i-1)*DTHETA+STHETA+PI)+TANK_C2[1];

      p.tag = ps.size()+1;
      p.r = xx,yy;
      p.dens = DENS;
      p.mass = PSEP*PSEP*DENS;
      p.h = H;
      p.v = 0.0,0.0;
      p.iam = boundary;
      p.norm1[0] = -cos((i-1)*DTHETA+STHETA+PI);
      p.norm1[1] = -sin((i-1)*DTHETA+STHETA+PI);
      if (i==NTANKPERIM) {
         p.concave = 0;
         p.norm2[0] = -cos(STHETA);
         p.norm2[1] = -sin(STHETA);
      }
      p.alpha = ALPHA;
      ps.push_back(p);

   }


      
   cout << "Total number of particles = " << ps.size() << endl;

   CglobalVars globals;


   vector<vector<double> > vprocDomain(globals.mpiSize);
   vector<Array<int,NDIM> > vprocNeighbrs(globals.mpiSize);
   vector<particleContainer > vps;
   vectInt split;
   split = 1,1;
   particleContainer pps;
   for (int i=0;i<ps.size();i++) {
      pps.push_back(ps[i]);
   }
   Nmisc::splitDomain(pps,split,vps,vprocDomain,vprocNeighbrs);



   //CdataLL *data = new CdataLL(ps,globals);
   //CsphIncompress sph(data);
   //CcustomOutput customOutput(data);
   //CcustomSim customSim(data,globals.time);


   cout << "Opening files for writing..."<<endl;
   Cio_data_vtk ioFile(filename.c_str(),&globals);
   cout << "Calculating Output stuff.."<<endl;
   //sph.calcOutputVars();
   //customOutput.calcOutput(0,&customSim,&ioFile);
   cout << "Writing Restart data to file..."<<endl;
   int nProc = product(split);
   for (int i=0;i<nProc;i++) {
      globals.mpiRank = i;
      ioFile.setFilename(filename.c_str(),&globals);
      ioFile.writeGlobals(0,&globals);
      ioFile.writeRestart(0,vps[i],&globals);
      globals.mpiRank = 0;
   }
   cout << "Writing Global data to file..."<<endl;
   //ioFile.writeGlobals(0,&globals);
   ioFile.writeDomain(0,vprocDomain,vprocNeighbrs);


   //Write restart file for John
   /*ofstream fo("restartJohn.dat");
   fo <<"h = "<<H<<" mass = "<<ps[0].mass<<" dens = "<<ps[0].dens<<" psep = "<<PSEP<<" Re = "<<REYNOLDS_NUMBER<<" kinematic viscosity = "<<VISCOSITY<<" n = "<<ps.size()<<endl;
   fo <<"r_x r_y v_x v_y flag(0=fluid,1=boundary)"<<endl;
   for (int i=0;i<ps.size();i++) {
      fo << ps[i].r[0]<<' '<<ps[i].r[1]<<' '<<ps[i].v[0]<<' '<<ps[i].v[1]<<' '<<ps[i].iam<<endl;
   }*/

   }

