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


#ifdef VAR_H_CORRECTION
      static double readOmega(Cparticle &p) {
         return p.omega;
      }
      static void writeOmega(Cparticle &p,double omega) {
         p.omega = omega;
      }
#endif

#ifdef VAR_H_CORRECTION2
      static vect readGradH(Cparticle &p) {
         return p.gradH;
      }
      static void writeGradH(Cparticle &p,vect gradH) {
         p.gradH = gradH;
      }
      static void sumGradH(Cparticle &p,vect gradH) {
         p.gradH += gradH;
      }
#endif

#ifdef LIQ_DEM
      static double readShepSum(Cparticle &p) {
         return p.shepSum;
      }
      static void writeShepSum(Cparticle &p,double shepSum) {
         p.shepSum = shepSum;
      }
      static void sumShepSum(Cparticle &p,double shepSum) {
         p.shepSum += shepSum;
      }
      static double readPorosity(Cparticle &p) {
         return p.porosity;
      }
      static double read_dPorositydt(Cparticle &p) {
         return p.dporositydt;
      }
      static void writePorosity(Cparticle &p,double porosity) {
         p.porosity = porosity;
      }
      static void write_dPorositydt(Cparticle &p,double dporositydt) {
         p.dporositydt = dporositydt;
      }
      static void sumPorosity(Cparticle &p,double porosity) {
         p.porosity += porosity;
      }
      static void sum_dPorositydt(Cparticle &p,double dporositydt) {
         p.dporositydt += dporositydt;
      }
      static vect readFdrag(Cparticle &p) {
         return p.fdrag;
      }
      static void writeFdrag(Cparticle &p,vect fdrag) {
         p.fdrag = fdrag;
      }
      static void sumFdrag(Cparticle &p,vect fdrag) {
         p.fdrag += fdrag;
      }
#endif


      static void calcAvForce(Cparticle &p,CglobalVars &g,vect &avForce) {
         avForce += p.f;
         g.tmp++;
      }

      static void setAvForce(Cparticle &p,CglobalVars &g,vect &avForce) {
         p.f = avForce;
      }

      /*
      static double readhErr(Cparticle &p) {
         return p.hErr;
      }
      static void writehErr(Cparticle &p,double hErr) {
         p.hErr = hErr;
      }
      */ 
      static vect readFp(Cparticle &p) {
         return p.fp;
      }
      static void writeFp(Cparticle &p,vect fp) {
         p.fp = fp;
      }
      static vect readFv(Cparticle &p) {
         return p.fv;
      }
      static void writeFv(Cparticle &p,vect fv) {
         p.fv = fv;
      }
      static vect readFb(Cparticle &p) {
         return p.fb;
      }
      static void writeFb(Cparticle &p,vect fb) {
         p.fb = fb;
      }
      static double readTmp(Cparticle &p) {
         return p.tmp;
      }
      static void writeTmp(Cparticle &p,double tmp) {
         p.tmp = tmp;
      }
#ifdef SLK
      static vectTensor readG(Cparticle &p) {
         return p.G;
      }
      static void writeG(Cparticle &p,vectTensor G) {
         p.G = G;
      }
#endif
#ifdef MY_VAR_RES
      static vect readGradV0(Cparticle &p) {
         return p.gradV[0];
      } 
      static void writeGradV0(Cparticle &p,vect a) {
         p.gradV[0] = a;
      }
      static vect readGradV1(Cparticle &p) {
         return p.gradV[1];
      } 
      static void writeGradV1(Cparticle &p,vect a) {
         p.gradV[1] = a;
      }
      static vect readGradV2(Cparticle &p) {
         return p.gradV[2];
      } 
      static void writeGradV2(Cparticle &p,vect a) {
         p.gradV[2] = a;
      }
#endif

#ifdef SPH_SMOOTHING
      static vect readBasis(Cparticle &p) {
         return p.basis;
      }
      static void writeBasis(Cparticle &p,vect basis) {
         p.basis = basis;
      }
#endif
      static vect readF(Cparticle &p) {
         return p.f;
      }
      static void writeF(Cparticle &p,vect f) {
         p.f = f;
      }
      static vect readVhat(Cparticle &p) {
         return p.vhat;
      }
      static void writeVhat(Cparticle &p,vect vhat) {
         p.vhat = vhat;
      }
      static vect readV(Cparticle &p) {
         return p.v;
      }
      static void writeV(Cparticle &p,vect v) {
         p.v = v;
      }
 
      static vect readR(Cparticle &p) {
         return p.r;
      }
      static void writeR(Cparticle &p,vect r) {
         p.r = r;
      }
      static double readDens(Cparticle &p) {
         return p.dens;
      }
      static void writeDens(Cparticle &p,double dens) {
         p.dens = dens;
         if (p.iam==sph) {
         for (int i=0;i<NDIM;i++) {
            if (PERIODIC[i]) {
               if (p.r[i]<RMIN[i]) {
                  p.dens += DENS_DROP[i];
               } else if (p.r[i]>RMAX[i]) {
                  p.dens -= DENS_DROP[i];
               }
            }
         }
         }
      }
      static double readH(Cparticle &p) {
         return p.h;
      }
      static void writeH(Cparticle &p,double h) {
         p.h = h;
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

#ifdef SMOOTHED_VISC_VELOCITY
      static vect readViscV(Cparticle &p) {
         return p.viscv;
      }
      static void writeViscV(Cparticle &p,vect viscv) {
         p.viscv = viscv;
      }
      static void smoothViscV(Cparticle &p, vector<Cparticle *> &neighbrs, CglobalVars &g) {
         if (p.iam == sphBoundary) {
            p.viscv = p.vhat;
            return;
         }
#ifdef SMOOTHED_VISC_VELOCITY_MLS
         double b0;
         vect bRest;
         vector<double> vWab;
         calcB_MLS(p,neighbrs,b0,bRest,vWab);
         int n = neighbrs.size();
         p.viscv = 0;
         for (int i=0;i<n;i++) {
            Cparticle *pn = neighbrs[i];
            p.viscv += (pn->vhat*pn->mass/pn->dens)*W_MLS(p.r-pn->r,vWab[i],b0,bRest);
         }
#endif
#ifdef SMOOTHED_VISC_VELOCITY_HAT
         p.viscv = 0;
         for (vector<Cparticle *>::iterator pNeighbr = neighbrs.begin();pNeighbr!=neighbrs.end();pNeighbr++) {
            Cparticle *pn = *pNeighbr;
            const double hav = 0.5*(p.h+pn->h);
            const double q = len(p.r-pn->r)/hav;
            
            double vK = pn->mass*W(q,hav)*2.0/(p.dens+pn->dens);
            p.viscv += pn->mass*W(q,hav)*2.0/(p.dens+pn->dens)*(pn->vhat-p.vhat);
         }
         p.viscv = p.vhat + EPSILON*p.viscv;
#endif
      }
#endif
      
#ifdef DIRECT_SMOOTHING
      static void driftVhatAndV(Cparticle &p,CglobalVars &g) {
         double dt = g.dt/2;
         p.vhat = 2*p.vhat - p.vhat0;
         p.v += p.f*dt;
      }

      static void driftAndSaveR(Cparticle &p,CglobalVars &g) {
         p.norm1 = p.r;
         driftR(p,g);
      }
      
      static void smoothVhatAndF(Cparticle &p, vector<Cparticle *> &neighbrs, CglobalVars &g) {
         p.vhat = 0;
         p.fhat = 0;  
         for (vector<Cparticle *>::iterator pNeighbr = neighbrs.begin();pNeighbr!=neighbrs.end();pNeighbr++) {
            Cparticle *pn = *pNeighbr;
            const double hav = 0.5*(p.h+pn->h);
            const double q = len(p.r-pn->r)/hav;
            
            double vK = pn->mass*W(q,hav)*2.0/(p.dens+pn->dens);
            p.vhat += vK*(pn->v-p.v);
            p.fhat += vK*(pn->f-p.f);
         }
         p.vhat = p.v + EPSILON*p.vhat;
         p.fhat = p.f + EPSILON*p.fhat;
      }
      
      static void calcVhatFromV(Cparticle &p, vector<Cparticle *> &neighbrs, CglobalVars &g) {
         vect oldVhat = p.vhat;
         p.vhat = 0;
         for (vector<Cparticle *>::iterator pNeighbr = neighbrs.begin();pNeighbr!=neighbrs.end();pNeighbr++) {
            Cparticle *pn = *pNeighbr;
            const double hav = 0.5*(p.h+pn->h);
            const double q = len(p.r-pn->r)/hav;
            
            double vK = pn->mass*W(q,hav)*2.0/(p.dens+pn->dens);
            p.vhat += vK*(pn->v-p.v);
         }
         p.vhat = p.v + EPSILON*p.vhat;
         g.residual += abs(oldVhat-p.vhat);
         g.DSTmpInt++;
      }

#endif           
#ifdef SPH_SMOOTHING
/*
      static void driftVhat(Cparticle &p,CglobalVars &g) {
         if (p.tag==1000) cout <<"p.vhat0 = "<<p.vhat0<<" p.vhat(1/2) = "<<p.vhat;
         p.vhat = 2.0*p.vhat - p.vhat0;
         if (p.tag==1000) cout <<" p.vhat(1) = "<<p.vhat<<endl;
      }
*/


/*
      static void driftRAndincrH(Cparticle &p,CglobalVars &g) {
         if (p.iam==sph) {
            p.norm = p.r;
            double dt = g.dt/2;
            p.r += dt*p.vhat;
         }
         
         p.v0[0] = p.h;
         p.h *= SMOOTH_HFAC_MULT;
      }

      static void updateVhat(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g) {
         if (neighbrs.size()>=NUM_SMOOTH_NEIGHBOURS) {
            cerr << "Error (updateVhat): too many neighbours. NUMNEIGHBOURS = "<<NUM_SMOOTH_NEIGHBOURS<<" neighbrs.size() = "<<neighbrs.size()<<" !!"<<endl;
            exit(-1);
         }
         vect sumVhat = 0.0;
         double sumW = 0;
         double tmpVal,hav;
         for (vector<Cparticle *>::iterator pNeighbr = neighbrs.begin();pNeighbr!=neighbrs.end();pNeighbr++) {
            Cparticle *pn = *pNeighbr;
            tmpVal = pn->mass*W(len(p.r-pn->r)/H,H)/DENS;
            sumVhat += tmpVal*pn->vhat;
            sumW += tmpVal;
         }
         p.vhat = (p.v + EPSILON*sumVhat)/(1.0 + EPSILON*sumW);
      }

      static void resetRandH(Cparticle &p,CglobalVars &g) {
         if (p.iam==sph) {
            double dt = g.dt/2;
            p.r = p.norm + dt*p.vhat;
         }
         p.h = p.v0[0];
      }
         
      
         
      static void resetH(Cparticle &p,CglobalVars &g) {
         p.h = p.v0[0];
      }

*/
      static void calcResidualsImplicitSmoothing(Cparticle &p,CglobalVars &g) {
         p.vhat += g.alpha*p.basis;
         p.residual -= g.alpha*p.pTimesA;
         g.residual2 += p.residual*p.residual/p.matA[0];  //using diagonal pre-conditioner
         //g.residual2 += p.residual*p.residual;  
      }

      static void updateVhatImplicitSmoothing(Cparticle &p,CglobalVars &g) {
         p.vhat += g.alpha*p.basis;
      }

      static void calcNewBasisImplicitSmoothing(Cparticle &p,CglobalVars &g) {
         p.basis = p.residual/p.matA[0] + g.beta*p.basis; //using diagonal pre-conditioner
         //p.basis = p.residual + g.beta*p.basis; 

      }

      static void confirmResidual2(Cparticle &p, vector<Cparticle *> &neighbrs, CglobalVars &g) {

         int n = neighbrs.size();
         vect residualSum = 0.5*(p.v + p.v0) - p.vhat;
         double hav;
         for (int i=0;i<n;i++) {
            Cparticle *pn = neighbrs[i];
            hav = 0.5*SMOOTH_HFAC_MULT*(p.h+pn->h);
            residualSum += EPSILON*pn->mass*2.0*(pn->vhat-p.vhat)*W(len(p.r-pn->r)/hav,hav)/(p.dens+pn->dens);
            //residualSum += EPSILON*pn->mass*(pn->vhat-p.vhat)*W(len(p.r-pn->r)/H,H)/DENS;
         }
         g.residual2 += residualSum*residualSum;
      }

      static void initImplicitSmoothing(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g) {
         if (neighbrs.size()>=NUM_SMOOTH_NEIGHBOURS) {
            cerr << "Error (initImplicitSmoothing): too many neighbours. NUMNEIGHBOURS = "<<NUM_SMOOTH_NEIGHBOURS<<" neighbrs.size() = "<<neighbrs.size()<<" !!"<<endl;
            exit(-1);
         }
         p.residual = 0.5*(p.v + p.v0);
         double sumForA = 0.0;
         double compA,hav;
         vector<double>::iterator pmatA = p.matA.begin();
         for (vector<Cparticle *>::iterator pNeighbr = neighbrs.begin();pNeighbr!=neighbrs.end();pNeighbr++) {
            pmatA++;
            Cparticle *pn = *pNeighbr;
            hav = 0.5*SMOOTH_HFAC_MULT*(p.h+pn->h);
            compA = EPSILON*pn->mass*2.0*W(len(p.r-pn->r)/hav,hav)/(p.dens+pn->dens);
            //compA = EPSILON*pn->mass*W(len(p.r-pn->r)/H,H)/DENS;
            sumForA += compA;
            *pmatA = -compA;
            p.residual += compA*pn->vhat;
         }
         p.matA[0] = 1.0+sumForA;
         p.residual -= p.matA[0]*p.vhat;
         p.basis = p.residual/p.matA[0];  //using diagonal pre-conditioner
         //p.basis = p.residual;  
         g.residual2 += p.residual*p.basis; 
      }

      static void calcAlphaDenomImplicitSmoothing(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g) {
         vector<double>::iterator pmatA = p.matA.begin();
         p.pTimesA = p.basis*(*pmatA);
         for (vector<Cparticle *>::iterator pNeighbr = neighbrs.begin();pNeighbr!=neighbrs.end();pNeighbr++) {
            pmatA++;
            Cparticle *pn = *pNeighbr;
            p.pTimesA += pn->basis*(*pmatA);
         }
         g.alphaDenom += p.pTimesA*p.basis;
      }

#endif            
      static void resetR(Cparticle &p,CglobalVars &g) {
         double dt = g.dt/2;
         p.r = p.norm1 + dt*p.vhat;
      }

#ifdef GRID_SMOOTHING
      static void setVtoVhat(Cparticle &p,CglobalVars &g) {
         p.v = p.vhat;
      }
      static void setVhatToV(Cparticle &p,CglobalVars &g) {
         p.vhat = p.v;
      }
      static void initGridParticles(Cparticle &p,CglobalVars &g) {
         double dt = g.dt/2;
         p.norm1 = p.r;
         p.r += dt*p.vhat;
      } 

      static void linearSmoothingP_TO_G(Cparticle &p,Array<CsmVertex*,NDIM> &vertices,CglobalVars &g) {
      //static void linearSmoothingP_TO_G(Cparticle &p,vector<CsmVertex *> &vertices,CglobalVars &g) {
         const vectInt coordsBL = 0;
         const vectInt coordsTR = 1;
         const vect vbase = (p.mass/(pow(SMOOTH_GRID_SIZE,NDIM)*p.dens))*p.v;
         //const vect vbase = (1.0/(pow(SMOOTH_GRID_SIZE,NDIM)))*p.v;
         const vect pos = (p.r-vertices(coordsBL)->pos)/(vertices(coordsTR)->pos-vertices(coordsBL)->pos); 
         for (Array<CsmVertex*,NDIM>::iterator i=vertices.begin();i!=vertices.end();i++) {
         //for (vector<CsmVertex*,NDIM>::iterator i=vertices.begin();i!=vertices.end();i++) {
            bool onBoundary = false;
            for (int j=0;j<NDIM;j++) {
               if (((*i)->pos[j]==RMIN[j]) || (abs((*i)->pos[j]-RMAX[j])<SMOOTH_GRID_SIZE*0.0001)) {
                  //cout <<"on boundary: i->pos = "<<i->pos<<endl;
                  onBoundary = true;
                  break;
               }
            }
            if (!onBoundary) {
               //double q = len(p.r-(*i)->pos)*2/SMOOTH_GRID_SIZE;
               vect v = vbase;
               vectInt coords = i.position();
               for (int j=0;j<NDIM;j++) {
                  if (coords[j]==0) {
                     v *= (1-pos[j]);
                  } else {
                     v *= pos[j];
                  }
               }
               (*i)->v += v;
               //i->v = 1;
               //(*i)->v += W(q,SMOOTH_GRID_SIZE/2)*p.mass*p.v/p.dens;
               (*i)->num++;
            }
         }
      }
   
      static void initGrid(CsmVertex &vertex,CglobalVars &g) {
         vertex.v = 0;
         vertex.num = 0;
      }
   
      static void finaliseGrid(CsmVertex &vertex,CglobalVars &g) {
         //if (vertex.num!=0) {
         //   vertex.v /= vertex.num;
            //cout <<"v = "<<vertex.v<<" num = "<<vertex.num<<endl;
         //}
         //vertex.v /= 4.0*pow(SMOOTH_GRID_SIZE,2);
      }

      //static void linearSmoothingG_TO_P(Cparticle &p,vector<CsmVertex *> &vertices,CglobalVars &g) {
      static void linearSmoothingG_TO_P(Cparticle &p,Array<CsmVertex*,NDIM> &vertices,CglobalVars &g) {
         const vectInt coordsBL = 0;
         const vectInt coordsTR = 1;
         vect pos = (p.r-vertices(coordsBL)->pos)/(vertices(coordsTR)->pos-vertices(coordsBL)->pos); 
         vect oldVhat = p.vhat;
         p.vhat = 0;
         //cout <<" BL = "<<vertices(coordsBL).pos<<" TR = "<<vertices(coordsTR).pos<<" pos = "<<pos;
         for (Array<CsmVertex*,NDIM>::iterator i=vertices.begin();i!=vertices.end();i++) {
            vectInt coords = i.position();
            vect v = (*i)->v;
            for (int j=0;j<NDIM;j++) {
               if (coords[j]==0) {
                  v *= (1-pos[j]);
               } else {
                  v *= pos[j];
               }
            }
            //cout <<" i->v = "<<i->v;
            p.vhat += v;
         }
         //cout <<" p.vhat = "<<p.vhat<<endl;
         g.itError2 += len2(oldVhat-p.vhat);
      }
#endif
      static void driftR(Cparticle &p,CglobalVars &g) {
         //if (p.tag==1000) cout << "driftR: dt= "<<g.dt/2<<" p.v = "<<p.v<<" p.f = "<<p.f<<endl;
         const double dt = g.dt/2;
#ifdef CHECK_FOR_NAN
         vect oldr = p.r;
#endif
#ifdef SLK
         p.r += p.dr[g.slkIndex0];
         p.thisRInc = p.dr[g.slkIndex0];
         const vect rInc = dt*p.vhat;
         p.dr[g.slkIndex0] = rInc;
         p.currR += rInc; 
#else
         p.r += dt*p.vhat;
#endif
#ifdef CHECK_FOR_NAN
         for (int i=0;i<NDIM;i++) {
            if (p.r[i]!=p.r[i]) 
               cout <<"ERROR: found nan for oldr = "<<oldr<<" tag = "<<p.tag<<endl;
         }
#endif
      }
      
      static void driftRAndKick(Cparticle &p,CglobalVars &g) {
         //if (p.tag==1000) cout << "driftAndKick: dt= "<<g.dt/2<<" p.v = "<<p.v<<" p.f = "<<p.f<<endl;
         const double dt = g.dt/2;

#ifdef SLK
         p.r += p.dr[g.slkIndex0];
         p.thisRInc = p.dr[g.slkIndex0];
         const vect rInc = dt*p.vhat;
         p.dr[g.slkIndex0] = rInc;
         p.currR += rInc; 
#else
         p.r += dt*p.vhat;
#endif
         if (p.iam!=sphBoundary) {
            p.v0 = p.v;
            p.vhat0 = p.vhat;
            p.v += dt*p.f;
            p.vhat += dt*p.f;
#ifdef FIXED_DEM
            if (p.iam==dem) {
               p.v = p.v0;
               p.vhat = p.vhat0;
            }
#endif
         }
         //if (g.sphStep < DAMP_STEPS) {
         //   p.v *= 0.99;
         //   p.vhat *= 0.99;
         //}
#ifdef SAVE_VHALF
         p.vHalf = p.v;
         p.rHalf = p.r;
#endif
      }

      static void driftRest(Cparticle &p,CglobalVars &g) {
         //if (p.tag==1000) cout << "driftRest: dt= "<<g.dt/2<<" p.dens = "<<p.dens<<" p.dddt = "<<p.dddt<<endl;
         const double dt = g.dt/2;
         const double pin = p.dens;
#ifdef SLK
         //const double massin = p.mass;
         //p.mass += p.mass*(p.dddt*dt - p.dRhoKernel)/p.dens;
         const double thisDrho = dt*p.dddt;
         //p.dens += p.dRhoKernel + dt*p.dMassDiff - p.drho[g.slkIndex0] + thisDrho;
         //p.dens += p.dRhoKernel + p.drho[g.slkIndex0] + thisDrho;
         p.dens += p.dRhoKernel;
         p.drho[g.slkIndex0] = thisDrho; 
         //p.dens += dt*p.dddt;
         //const double dm = dt*p.dMassDiff;
         //p.dens += p.dens*dm/p.mass;
         //p.mass += dm;
#else
         p.dens += dt*p.dddt;
#endif
         //p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
#ifndef CONST_H
         p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
#endif
         //p.u += dt*p.dudt;
         //p.h += dt*p.dhdt;
         //p.alpha += dt*p.dalphdt;
      }

      static void kick(Cparticle &p,CglobalVars &g,double halfDt0,double halfDt1) {
         //if (p.tag==1000) cout << "kick: dt= "<<dt<<" p.v0 = "<<p.v0<<" p.f= "<<p.f<<" p.fp = "<<p.fp<<" p.fv = "<<p.fv<<endl;
         double dt = halfDt0+halfDt1;
#ifdef CHECK_FOR_NAN
         vect oldv = p.v;
#endif
         if ((p.iam == sph)||(p.iam==dem)||(p.iam==immersedDem)) {
            p.v = p.v0 + dt*p.f;
            p.vhat = p.vhat0 + dt*p.f;
            if (g.time < DAMPTIME) {
               p.v *= 0.98;
               p.vhat *= 0.98;
            }
#ifdef FIXED_DEM
            if (p.iam==dem) {
               p.v = p.v0;
               p.vhat = p.vhat0;
            }
#endif
#ifdef AVE_VELOCITY
            p.aveV = 0;
            for (int i=0;i<AVE_VELOCITY_N-1;i++) {
               p.pastVs[i] = p.pastVs[i+1];
               p.aveV += p.pastVs[i];
            }
            p.pastVs[AVE_VELOCITY_N-1] = p.v;
            p.aveV += p.v;
            p.aveV /= AVE_VELOCITY_N;
#endif
         } 
#ifdef CHECK_FOR_NAN
         for (int i=0;i<NDIM;i++) {
            if (p.v[i]!=p.v[i]) 
               cout <<"ERROR: found nan for oldv = "<<oldv<<" tag = "<<p.tag<<endl;
         }
#endif
      }

      static void initSumsMiddle(Cparticle &p,CglobalVars &g) {
         //if (p.tag==1000) cout << "initSums:"<<endl;
         p.f = 0.0;
         p.fp = 0.0;
         p.fv = 0.0;
         p.fb = 0.0;
         p.ff = 0.0;
         p.dudt = 0.0;
         p.deViscFdt = 0.0;
         p.deViscBdt = 0.0;
         p.deBForcedt = 0.0;
         p.deFFdt = 0.0;
         //p.dhdt = 0.0;
         //p.dalphdt = 0.0;
         //p.colour = 0;
#ifdef SLK
         p.dRhoKernel = 0;
         p.dMassDiff = 0;
#endif
      }

      static void initSumsEnd(Cparticle &p,CglobalVars &g) {
         p.dddt = 0.0;
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

#ifdef VAR_H_CORRECTION
         const vect gradWa = gradW(pa,pb,dr,r/pa.h,pa.h);
#else
         const vect gradWa = 0.5*(gradW(pa,pb,dr,r/pa.h,pa.h)+gradW(pa,pb,dr,r/pb.h,pb.h));
#endif
#ifdef SLK
         const vect drInc = pa.thisRInc-pb.thisRInc;
         //if ((pb.iam==sphBoundary)&&((len2(pb.thisRInc)!=0)||any(pb.currR!=pb.r))) {
         //   cout << "CRAP: thisRInc = "<<pb.thisRInc<<", currR-r = "<<pb.currR-pb.r<<endl;
         //   exit(-1);
         //}
         pa.dRhoKernel += pb.mass*dot(drInc,gradWa);

         vect transGradWa = product(0.5*(pa.G+pb.G),gradWa);
         
         double newDr = 0.0;
         vect newDx = pa.currR-pb.currR;
         double newVdr = dot(newDx,dv);
         double newR = len(newDx);
         if (newR!=0.0) newDr = 1/newR;
         double viss = newVdr*newDr;
         double vsig = 2.0*SPSOUND + 2.0*abs(viss);
         const vect massDiff = 2.0*ALPHA*vsig*newDx*newDr*(pa.mass-pb.mass)/(pa.dens+pb.dens);
         //double massDiff = 2.0*ALPHA*newR*vsig*(pa.dens-pb.dens)/(pa.dens+pb.dens);
         pa.dMassDiff += pb.mass*dot(massDiff,transGradWa);
         
         const double vdr = dot(dv,transGradWa);
#else
         const double vdr = dot(dv,gradWa);
#endif

#ifdef DDENS_VARIANT
         const double dddtInc = pb.mass*vdr*(pa.dens/pb.dens)*sqrt(pb.press/pa.press);
#else
         const double dddtInc = pb.mass*vdr;
#endif

         
#ifdef LIQ_DEM_DDDT_VER2
         pa.dddt += (pa.porosity/pb.porosity)*(dddtInc + pb.mass*W(r/pa.h,pa.h)*(pa.dporositydt/pa.porosity - pb.dporositydt/pb.porosity));
#else
         pa.dddt += dddtInc;
#endif
         //cout <<"pa.dddt = "<<dddtInc<<"pb.mass = "<<pb.mass<<" vdr "<<vdr<<" q = "<<q<<"hav = "<<hav<<"Fa = "<<Fa<<" pa.r = "<<pa.r<<" pb.r "<<pb.r<<endl;

#ifdef DENS_DIFFUSE
         double rr = 0.0;
         if (r!=0.0) rr = 1/r;
         double viss = vdr*rr;
         //cout << "abs(viss) = "<<abs(viss)<<" avspsound = "<<0.5*(pa.spsound+pb.spsound)<<endl;
         double vsig = 2.0*SPSOUND + 2.0*abs(viss);
         //double vsig = 2.0*SPSOUND;
         //if (vsig > pa.maxvsig) pa.maxvsig = vsig;
         //double densDiff = 2.0*ALPHA*r*vsig*(pa.dens-pb.dens)/(pa.dens+pb.dens);
         double densDiff = 4.0*ALPHA*vsig*(pa.dens-pb.dens)/(pa.dens+pb.dens);
         pa.dddt += pb.mass*densDiff*dot(dr*rr,gradWa);
#endif

         //if (pa.tag==1000) cout << "calcDddtDudt: dr = "<<dr<<" dv = "<<dv<<" r/hav = "<<r/hav<<" hav = "<<hav<<" vdr = "<<vdr<<" Fa = "<<Fa<<" pb.mass = "<<pb.mass<<" pa.dddt = "<<pa.dddt<<" W = "<<W(r/hav,hav)<<endl;
         //pa.dudt += pa.pdr2*dddtInc;
         //pa.colour += 1;
      }

//void CsphIncompress:init_particles(Cparticle &p) {
//}
   
#ifdef LIQ_DEM
      static void correctDddt(Cparticle &p,CglobalVars &g) {
         p.dddt -= p.dens*dot(p.vhat,p.gradPorosity);
         //p.dddt *= 1.0/p.porosity;
      }
#endif

      static void calcPressSpsoundPdr2(Cparticle &p,CglobalVars &g) {
#ifdef INCL_THERM_ENERGY
         p.press = PRB*(1.0+p.u)*(pow(p.dens/REFD,7) - 1.0);
#else
#ifdef LIQ_DEM
#ifdef LIQ_DEM_DENS_NORMAL
         p.press = PRB*(pow(p.dens/REFD,7) - 1.0);
#else
         p.press = PRB*(pow(p.dens/(REFD*p.porosity),GAMMA) - 1.0);
#endif
         //p.press = PRB*(pow(p.dens/REFD,7) - 1.0);
#else
         p.press = PRB*(pow(p.dens/REFD,7) - 1.0);
         //p.press = RMAX[2]-p.r[1];
         //cout << "press = "<<p.press<<" PRB = "<<PRB<<" dens = "<<p.dens<<endl;
         //p.press = PRB*pow(p.dens/REFD,1);
#endif
#endif
         //p.press = PRB*(pow(p.dens/REFD,7) - 0.99);
         //p.press = pow(SPSOUND,2)*p.dens;
         p.pdr2 = p.press/(p.dens*p.dens);
#ifdef LIQ_DEM
#ifdef LIQ_DEM_SEPARATE_DRAGS
         p.pdr2 *= p.porosity;
#else
#ifdef LIQ_DEM_TEST
         p.pdr2 *= p.porosity;
#endif
#endif
#endif
#ifdef LIQ_DEM_DENS_NORMAL
         p.pdr2 /= p.porosity*p.porosity;
#endif
      }

      
      static void calcForce(Cparticle &pa, Cparticle &pb,CglobalVars &g) {

         vect dx = pa.r-pb.r;
         const double r = len(dx);
#ifdef SLK
         const double hav = 0.5*(pa.h+pb.h);
         const vect gradWa = product(0.5*(pa.G+pb.G),gradW(pa,pb,dx,r/hav,hav));
         const vect gradWb = gradWa;
#else
#ifdef VAR_H_CORRECTION
         const vect gradWa = gradW(pa,pb,dx,r/pa.h,pa.h);
         const vect gradWb = gradW(pa,pb,dx,r/pb.h,pb.h);
#else
         const vect gradWa = 0.5*(gradW(pa,pb,dx,r/pa.h,pa.h)+gradW(pa,pb,dx,r/pb.h,pb.h));
         const vect gradWb = gradWa;
#endif
#endif


         if ((pa.iam == sph)&&(pb.iam != dem)&&(pb.iam != demBoundary)) {
#ifdef MORRIS_SPH_BOUNDARY
            vect dv;
            if (pb.iam == sphBoundary) {
               vect norm;
               calcNormal(pa,pb,dx,r,norm);
               calcMorrisDv(pa,pb,norm,dx,dv);
            } else {
               dv = pa.v-pb.v;
            }
#else
#ifdef SMOOTHED_VISC_VELOCITY
            const vect dv = pa.viscv-pb.viscv;
#else
#ifdef DIRECT_SMOOTHING
            const vect dv = pa.vhat-pb.vhat;
#else
            //for boundary particle vhat is used for its movement, and v is used for the boundary velocity
            //const vect dv = pa.v-pb.v;    /*only use this for viscous force*/
            const vect dv = pa.vhat-pb.vhat;    /*only use this for viscous force*/
#endif
#endif
#endif

#ifdef SLIP_BOUNDARIES
            if (!ifBoundary(pb)) calcViscForce(pa,pb,g,0.5*(gradWa+gradWb),dv,r);
#else
            calcViscForce(pa,pb,g,0.5*(gradWa+gradWb),dv,r);
#endif
               

         }
        
#ifdef NO_ANTICLUMPING
         double kdwPowNeps = 0.0;
#else
         const double hav = 0.5*(pa.h+pb.h);
         const double q = r/hav;
         const double kdw = K(q,hav)/K(1.0/HFAC,1.0);
         const double kdwPowNeps = pow(kdw,NEPS);
#endif
         //cout << "wdwPowNeps = "<<wdwPowNeps<<"W1 = "<<W1<<"W(q) = "<<W(q,H)<<" q = "<<q<<endl;

         if (ifBoundary(pb)) {
#ifdef RADIAL_BOUNDARY
            calcRadialBoundaryForces2D(pa,pb,g);
#else
            vect normp,normt;
            calcNormal(pa,pb,dx,r,normp,normt);
            calcBoundaryForces(pa,pb,g,normp,normt,kdwPowNeps);
#endif
         } else if ((pb.iam == sph)||(pb.iam == sphBoundary)) {
            calcPressForce(pa,pb,g,gradWa,gradWb,kdwPowNeps);
         }
      }

   private:
      CdataLL *data;

};

#endif
