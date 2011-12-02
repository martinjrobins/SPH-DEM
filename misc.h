#ifndef MISC_H
#define MISC_H

#include <list>
#include <vector>
#include <blitz/array.h>
#include "dataLL.h"
#include "customConstants.h"
#include "vect.h"
#include "globalVars.h"
#include "particle.h"
#include "io_data_vtk.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <blitz/array.h>
#include <complex>

#include <ext/slist>
namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later

//using namespace std;

typedef complex<double> myComplex;
namespace Nmisc {
   
   class CpairData {
      public:
         int p1,p2;
         double initR2,initT,finalR2;
         ofstream *outFile;
   };


   inline double alphaToVisc(const double alpha,const double h,const double spsound) {
      return (15.0/112.0)*alpha*h*spsound;
   }
   inline double viscToAlpha(const double visc,const double h,const double spsound) {
      return visc*(112.0/15.0)/(spsound*h);
   }
   inline double viscToH(const double visc,const double alpha,const double spsound) {
      return visc*(112.0/15.0)/(spsound*alpha);
   }

   int coordsToNum(vectInt coords,vectInt split); 
   vectInt numToCoords(int num,vectInt split); 
   void setupEntireDomain(double procDomain[NDIM*2],Array<int,NDIM> &procNeighbrs);
   void splitDomain(particleContainer &ps,vectInt &split,vector<particleContainer > &vps,vector<vector<double> > &vprocDomain,vector<Array<int,NDIM> > &vprocNeighbrs);
   void divideParticlesToCPUs(particleContainer &ps,vector<particleContainer > &vps,vector<vector<double> > &vprocDomain);
   void calcLyap(CdataLL *d1,CdataLL *d2,Cio_data_vtk *io,double dt);
   void calcPoincare(CdataLL *d1,CdataLL *d2,Cio_data_vtk *io,int numParticles,int numPeriods);
   void pairDoubling(CdataLL *data,Sgi::slist<CpairData> &pairData,vector<ofstream *> &openFiles,gsl_rng *rng,Cio_data_vtk *io);

   template<bool ifFunct(Cparticle &p)>
   Array<myComplex,NDIM> *calcVelocitySpectrum(CdataLL *data, const int numModes);

   template<bool ifFunct(Cparticle &p)>
   void outputVelocitySpectrum(CdataLL *data,Cio_data_vtk *io);

   template<bool ifFunct(Cparticle &p),double specFunct(Cparticle &p,CglobalVars &g)>
   Array<myComplex,NDIM> *calcSpectrum(CdataLL *data, const int numModes);

   template<bool ifFunct(Cparticle &p),void specFunct(Cparticle &p,CglobalVars &g,vectInt &coords)>
   void setSpectrum(CdataLL *data, const int numModes);

   template<bool ifFunct(Cparticle &p),double specFunct(Cparticle &p,CglobalVars &g)>
   void outputSpectrum(CdataLL *data,Cio_data_vtk *io,string name);

   template<bool ifFunct(Cparticle &p)>
   void getRandomPair(particleContainer &ps,gsl_rng *rng,Cparticle **p1,Cparticle **p2);

   template<bool ifFunct(Cparticle &p)>
   void getRandomNeighbrPair(CdataLL *data,gsl_rng *rng,Cparticle **p1,Cparticle **p2);

   template<bool ifFunct(Cparticle &p)>
   void calcVelocityStructure(CdataLL *data,gsl_rng *rng,Cio_data_vtk *io);

   template<bool ifFunct(Cparticle &p)>
   void calcFullVelocityStructure(CdataLL *data,Cio_data_vtk *io);
    
   template<bool ifFunct(Cparticle &p)>
   void calcPositionStructure(CdataLL *data,gsl_rng *rng,Cio_data_vtk *io);


   void boundaryLine(vector<Cparticle> &ps,vect *points,vect *normals,bool *concave,int numPoints,enum iamTypes iam);

   void boundaryPlane(vector<Cparticle> &ps,vect side1,vect side2,vect neighbrNorm1, bool smoothEdge1, vect neighbrNorm2, bool smoothEdge2,vect neighbrNorm3, bool smoothEdge3,vect neighbrNorm4, bool smoothEdge4,vect initPos,double thisPsep,bool flipNorm);
   void boundaryCircle(vector<Cparticle> &ps,const vect &origin,const double radiusMin,const double radiusMax,const vect normal);
   void boundaryCylinderNoTopBottom(vector<Cparticle> &ps,const vect &origin,const double radius,const double height);

   void boundaryCylinder(vector<Cparticle> &ps,const vect &origin,const double radius,const double height);
   void sphCylinder(vector<Cparticle> &ps,const vect &origin);


   void addGridDEMparticles(vector<Cparticle> &ps,const double minRange[NDIM],const double maxRange[NDIM],const double porosity);
   void addRandomDEMparticles(CdataLL *data,const double minRange[NDIM],const double maxRange[NDIM],const double porosity);
   void addRandomDEMparticlesCyl(CdataLL *data,const vect origin,const double radius,const double height,int n);
}

#endif
