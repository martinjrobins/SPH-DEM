#ifndef SPHINCOMPRESS_H
#define SPHINCOMPRESS_H

#include <sys/time.h>
#include "dataLL.impl.h"
#include "vect.h"
#include "particle.h"
#include "sph.h"
#include "customConstants.h"
#include "globalVars.h"


using namespace Nsph;

class CsphIncompress {
   public:
      CsphIncompress(CdataLL *_data) {
         data = _data;
      };
      
      void start();
      void middle();
      void end();
      
      static void calcGridQuantities(Cparticle &gridP,Cparticle &p,CglobalVars &g) {
         double Wav = W(len(gridP.r-p.r)/p.h,p.h);
         double vol = p.mass/p.dens;
         gridP.dens += p.mass*Wav;
         gridP.u += vol*p.u*Wav;
         gridP.v += vol*p.v*Wav;
         gridP.vort += vol*p.vort*Wav;
      }

      void renderToGrid(vector<Cparticle> &ps,vectInt &gridDims) {
         data->functionOverGrid<calcGridQuantities>(ps,gridDims);
      }
      
      static void calcOutput(Cparticle &p,CglobalVars &g) {
         if (p.iam == sph) {
            g.linMom += p.v*p.mass;
            //TODO: this is 2D specific!
            g.angMom += (p.r[0]*p.v[1]-p.r[1]*p.v[0])*p.mass;
            g.eKE += 0.5*p.mass*len2(p.vhat);
            g.eViscF += p.mass*p.eViscF;
            g.edViscFdt += p.mass*p.deViscFdt;
            g.eViscB += p.mass*p.eViscB;
            g.eBForce += p.mass*p.eBForce;
            g.eFF += p.mass*p.eFF;
            g.edFFdt += p.mass*p.deFFdt;
            g.nSph++;
         }
         double v = len(p.v);
         if (v > g.maxV) g.maxV = v;
         double f = len(p.f);
         if (f > g.maxF) g.maxF = f;
         double ff = len(p.ff);
         if (ff > g.maxFF) g.maxFF = ff;
         g.rmsFF += len2(ff);
         g.aveDens += p.dens;
         g.eElast += p.mass*p.u;
         g.edElastdt += p.mass*p.dudt;
         double denfac = p.dens/REFD;
         denfac = PRB*((pow(denfac,6)-1)/6 + 1/denfac - 1);
         g.eElastExact += p.mass*denfac/REFD;
      }

      static void calcVarDens(Cparticle &p,CglobalVars &g) {
         g.varDens += pow(p.dens-g.aveDens,2);
      }

      static void calcAveDensFromMass(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g) {
         double oldDens = p.dens;
         calcDensity(p,neighbrs,g);
         g.aveDensFromMass += pow(p.dens-oldDens,2);
         p.tmp = p.dens;
         p.dens = oldDens;
      }

      void calcOutputVars() {
         data->globals.linMom = 0;
         data->globals.angMom = 0;
         data->globals.eKE = 0;
         data->globals.eFF = 0;
         data->globals.edFFdt = 0;
         data->globals.eViscF = 0;
         data->globals.edViscFdt = 0;
         data->globals.eViscB = 0;
         data->globals.eElast = 0;
         data->globals.edElastdt = 0;
         data->globals.eElastExact = 0;
         data->globals.eBForce = 0;
         data->globals.maxV = 0;
         data->globals.maxF = 0;
         data->globals.maxFF = 0;
         data->globals.aveDens = 0;
         data->globals.rmsFF = 0;
         data->globals.varDens = 0;
         data->globals.aveDensFromMass = 0;
         data->globals.nSph = 0;
         data->traverse<calcOutput,ifSphOrSphBoundary>();
         data->neighboursGroup<calcAveDensFromMass,ifSphOrSphBoundary>();
         data->globals.n = data->getParticles()->size();
         double tmpData[5];
         tmpData[0] = data->globals.n;
         tmpData[1] = data->globals.aveDens;
         tmpData[2] = data->globals.aveDensFromMass;
         tmpData[3] = data->globals.nSph;
         tmpData[4] = data->globals.rmsFF;
         data->sumOverProcs(tmpData,5);
         data->globals.n = (int)tmpData[0];
         data->globals.nSph = (int)tmpData[3];
         data->globals.aveDens = tmpData[1]/tmpData[0];
         data->globals.aveDensFromMass = sqrt(tmpData[2]/tmpData[0]);
         data->globals.rmsFF = sqrt(tmpData[4]/tmpData[0]);
         data->traverse<calcVarDens,ifSphOrSphBoundary>();
         data->sumOverProcs(&(data->globals.varDens),1);
         data->globals.varDens /= tmpData[0];
         data->globals.eTotal = data->globals.eKE 
                              + data->globals.eElast
                              + data->globals.eViscF 
                              + data->globals.eViscB 
                              + data->globals.eBForce
                              + data->globals.eFF;
      }

      static void driftR(Cparticle &p,CglobalVars &g) {
         const double dt = g.dt/2;
         p.r += dt*p.vhat;
      }
      
      static void driftRAndKick(Cparticle &p,CglobalVars &g) {
         const double dt = g.dt/2;
         p.r += dt*p.vhat;
         if (p.iam!=sphBoundary) {
            p.v0 = p.v;
            p.vhat0 = p.vhat;
            p.v += dt*p.f;
            p.vhat += dt*p.f;
         }         
      }

      static void driftRest(Cparticle &p,CglobalVars &g) {
         const double dt = g.dt/2;
         const double pin = p.dens;
         p.dens += dt*p.dddt;
#ifndef CONST_H
         p.h *= pow(pin/p.dens,1.0/NDIM);
#endif
      }

      static void kick(Cparticle &p,CglobalVars &g,double halfDt0,double halfDt1) {
         double dt = halfDt0+halfDt1;
         if ((p.iam == sph)||(p.iam == dem)) {
            p.v = p.v0 + dt*p.f;
            p.vhat = p.vhat0 + dt*p.f;

            if (g.time < DAMPTIME) {
               p.v *= 0.998;
               p.vhat *= 0.998;
            }
      }

      static void initSums(Cparticle &p,CglobalVars &g) {
         p.f = 0.0;
         p.fp = 0.0;
         p.fv = 0.0;
         p.fb = 0.0;
         p.ff = 0.0;
         p.dddt = 0.0;
         p.dudt = 0.0;
         p.deViscFdt = 0.0;
         p.deViscBdt = 0.0;
         p.deBForcedt = 0.0;
         p.deFFdt = 0.0;
      }

      static void calcEnergies(Cparticle &p,CglobalVars &g) {
         if (p.iam==sph) {
            p.eViscB = 0;
            p.eBForce = 0;
            p.deViscFdt = -dot(p.vhat,p.fv);
            p.eViscF += g.dt*p.deViscFdt;
            p.deFFdt = -dot(p.vhat,p.ff);
            p.eFF += g.dt*p.deFFdt;
         }
         calcPressSpsoundPdr2(p,g);
         p.dudt = p.pdr2*p.dddt;
         p.u += g.dt*p.dudt;
      }
      
      static void calcDddtDudt(Cparticle &pa, Cparticle &pb,CglobalVars &g) {
         if (ifBoundary(pb)||ifDemOrDemBoundary(pb)) return;
         vect dr = pa.r-pb.r;
         vect dv = pa.vhat-pb.vhat;
         const double r = len(dr);
         const double hav = 0.5*(pa.h+pb.h);
         const double q = r/hav;
         const vect gradWa = gradW(pa,pb,dr,q,hav);
         const double vdr = dot(dv,gradWa);
         const double dddtInc = pb.mass*vdr;
         pa.dddt += dddtInc;
      }

      static void calcPressSpsoundPdr2(Cparticle &p,CglobalVars &g) {
         p.press = PRB*(pow(p.dens/REFD,7) - 1.0);
         p.pdr2 = p.press/(p.dens*p.dens);
      }

      
      static void calcForce(Cparticle &pa, Cparticle &pb,CglobalVars &g) {
         vect dx = pa.r-pb.r;
         const double hav = 0.5*(pa.h+pb.h);
         const double r = len(dx);
         const double q = r/hav;
         const vect gradWa = gradW(pa,pb,dx,q,hav);

         if ((pa.iam == sph)&&(pb.iam != dem)&&(pb.iam != demBoundary)) {
            //for boundary particle vhat is used for its movement, and v is used for the boundary velocity
            //const vect dv = pa.v-pb.v;    /*only use this for viscous force*/
            const vect dv = pa.vhat-pb.vhat;    /*only use this for viscous force*/
            calcViscForce(pa,pb,g,gradWa,dv,r);
         }
        
#ifdef NO_ANTICLUMPING
         double kdwPowNeps = 0.0;
#else
         const double kdw = K(q,hav)/K(1.0/HFAC,1.0);
         const double kdwPowNeps = pow(kdw,NEPS);
#endif

         if (ifBoundary(pb)) {
            vect normp,normt;
            calcNormal(pa,pb,dx,r,normp,normt);
            calcBoundaryForces(pa,pb,g,normp,normt,kdwPowNeps);
         } else if ((pb.iam == sph)||(pb.iam == sphBoundary)) {
            calcPressForce(pa,pb,g,gradWa,kdwPowNeps);
         }
      }

   private:
      CdataLL *data;

};

#endif
