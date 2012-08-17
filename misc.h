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
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
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
   void boundaryCircle(vector<Cparticle> &ps,const vect &origin,const double radiusMin,const double radiusMax,const vect normal,bool inner_corner,bool outer_corner);
   void boundaryCylinderNoTopBottom(vector<Cparticle> &ps,const vect &origin,const double radius,const double height);
   void sphBoundaryCylinderNoTopBottom(vector<Cparticle> &ps,const vect &origin,const double radius,const double height);

   void boundaryCylinder(vector<Cparticle> &ps,const vect &origin,const double radius,const double height);
   void sphCylinder(vector<Cparticle> &ps,const vect &origin);

   /*
    * Add regular grid of dem particles within the rectangular block defined by
    * the points minRange and maxRange. The grid spacing is calculated using the 
    * desired porosity
    */
   void addGridDEMparticles(vector<Cparticle> &ps,const double minRange[NDIM],const double maxRange[NDIM],const double porosity);

   /*
    * Add dem particles randomly within the rectangle defined by the points
    * minRange and maxRange. The number of particles inserted is calculated 
    * using the desired porosity
    */
   void addRandomDEMparticles(CdataLL *data,const double minRange[NDIM],const double maxRange[NDIM],const double porosity);

   /*
    * Add dem particles randomly within a vertical cylinder. The origin point, radius
    * and height define the cylinder dimensions. The number of particles
    * inserted is calculated using the desired porosity
    */
   void addRandomDEMparticlesCyl(CdataLL *data,const vect origin,const double radius,const double height,int n);
}


#ifdef DEFORMATION_MATRIX
#include "sph.h"
namespace Nmisc {
inline void calcDeformationMatrix(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g,CdataLL &d2) {
   const int n = neighbrs.size();
   if (n<NDIM) {
      for (int j=0;j<NDIM;j++) {
         for (int i=0;i<n;i++) {
            p.defMat[i*NDIM+j] = 0;
         }
      }
      return;
   }
   gsl_matrix *X = gsl_matrix_alloc(n,NDIM);
   gsl_vector *y = gsl_vector_alloc(n);

   gsl_vector *c = gsl_vector_alloc(NDIM);
   gsl_matrix *cov = gsl_matrix_alloc(NDIM,NDIM);
            
   gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n,NDIM);

   vector<vect> dx2(n);
   //double shepSum = 0;
   for (int i=0;i<n;i++) {
      Cparticle *pn = neighbrs[i];
      vect dx = pn->r-p.r;
      /*
      const double r = len(dx);
      const double dvol = pn->mass/pn->dens;
      const double q = r/pn->h;
      const double Wab = Nsph::W(q,pn->h);
      const double dvWab = dvol*Wab;
      shepSum += dvWab;
      */
      Cparticle *pn2 = d2.getParticleFromTag(pn->tag);
      if (pn2==NULL) continue;
      for (int j=0;j<NDIM;j++) {
         if (PERIODIC[j]) {
            if (abs(pn2->r[j]-p.r[j]+RMAX[j]-RMIN[j])<abs(pn2->r[j]-p.r[j])) {
               (dx2[i])[j] = pn2->r[j]-p.r[j]+RMAX[j]-RMIN[j];
            } else if (abs(pn2->r[j]-p.r[j]-RMAX[j]+RMIN[j])<abs(pn2->r[j]-p.r[j])) {
               (dx2[i])[j] = pn2->r[j]-p.r[j]-RMAX[j]+RMIN[j];
            } else {
               (dx2[i])[j] = pn2->r[j]-p.r[j];
            }
         } else {
            (dx2[i])[j] = pn2->r[j]-p.r[j];
         }
         gsl_matrix_set(X,i,j,dx[j]);
      }
   }
   



   for (int j=0;j<NDIM;j++) {
      for (int i=0;i<n;i++) {
         gsl_vector_set(y,i,(dx2[i])[j]);
      }
            
      double chisq;
      gsl_multifit_linear(X,y,c,cov,&chisq,work);

      for (int i=0;i<NDIM;i++) {
         p.defMat[i*NDIM+j] = gsl_vector_get(c,i);
      }
   }
   gsl_multifit_linear_free(work);
   gsl_vector_free(c);
   gsl_vector_free(y);
   gsl_matrix_free(X);
   gsl_matrix_free(cov);
}



inline void calcRightCauchyGreenTensor(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g,CdataLL &d2) {

   calcDeformationMatrix(p,neighbrs,g,d2);

   // multiply transverse of coeff by itself to find Right Cauchy-Green
   // Deformation Tensor
   for (int i=0;i<NDIM;i++) {
      for (int j=0;j<NDIM;j++) {
         int index = i*NDIM+j;
         p.RCGMat[index] = 0;
         for (int k=0;k<NDIM;k++) {
            p.RCGMat[index] += p.defMat[k*NDIM+i]*p.defMat[k*NDIM+j];
         }
      }
   }
   double data[NDIM*NDIM];
   for (int i=0;i<NDIM;i++) {
      for (int j=0;j<NDIM;j++) {
         data[i*NDIM+j] = p.defMat[i*NDIM+j];
      }
   }

   //calculate eigenvalues of resultant matrix
   gsl_matrix_view m = gsl_matrix_view_array (data, NDIM, NDIM);
   gsl_vector *eval = gsl_vector_alloc (NDIM);
   gsl_matrix *evec = gsl_matrix_alloc (NDIM, NDIM);
   gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (NDIM);
   gsl_eigen_symmv(&m.matrix, eval, evec, w);
   gsl_eigen_symmv_free(w);
   gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);


   p.eval1 = gsl_vector_get(eval, NDIM-1);
   p.eval2 = gsl_vector_get(eval, NDIM-2);
   gsl_vector_view evec1 = gsl_matrix_column(evec,NDIM-1);
   for (int i=0;i<NDIM;i++) {
      p.evec1[i] = gsl_vector_get(&evec1.vector,i);
   }

   gsl_vector_free(eval);
   gsl_matrix_free(evec);

}

inline void calcParticleLyap(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g,CdataLL &d2) {
   calcRightCauchyGreenTensor(p,neighbrs,g,d2);
   double DLE = 1.0/(2.0*abs(d2.globals.time - g.time))*log(abs(p.eval1));
   p.tmp = DLE;
}

}
#endif



#endif
