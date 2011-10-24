
#ifndef VECT_H
#define VECT_H

#include "customConstants.h"
#include <blitz/tinyvec-et.h>
#include <gsl/gsl_linalg.h>
#include <blitz/array.h>


using namespace blitz;

typedef TinyVector<double,NDIM> vect;
typedef TinyVector<int,NDIM> vectInt;
typedef TinyVector<double,NDIM*NDIM> vectTensor;


#ifdef _2D_
inline vect cross(const vect &a,const vect &b) {
   //const TinyVector<double,NDIM> _a = a;
   //const TinyVector<double,NDIM> _b = b;

   //return cross(static_cast<TinyVector<double,NDIM> >(a),static_cast<TinyVector<double,NDIM> >(b));
   //return blitz::cross(_a,_b);
   return a[0]*b[1]-a[1]*b[0];
}
#endif

inline vectTensor outerProduct(vect a, vect b) {
   vectTensor returnT;
   for (int i=0;i<NDIM;i++) {
      for (int j=0;j<NDIM;j++) {
         returnT[i*NDIM+j] = a[i]*b[j];
      }
   }
   return returnT;
}

inline vect product(vectTensor a, vect b) {
   vect returnVect = 0.0;
   for (int i=0;i<NDIM;i++) {
      for (int j=0;j<NDIM;j++) {
         returnVect[i] += a[i*NDIM+j]*b[j];
      }
   }
   return returnVect;
}

inline void inverseSym(vectTensor &a_) {
#ifdef _2D_
   const double a = a_[0];
   const double b = a_[1];
   const double c = a_[2];
   const double d = a_[3];
   double deti =(a*d-b*c);
   if (deti != 0.0) {
      a_[0] = d;
      a_[1] = -b;
      a_[2] = -c;
      a_[3] = a;
      a_ = a_/deti;
   } else {
      a_[0] = 1.0;
      a_[1] = 0.0;
      a_[2] = 0.0;
      a_[3] = 1.0;
   }
#else
   const double a = a_[0];
   const double b = a_[1];
   const double c = a_[2];
   const double d = a_[3];
   const double e = a_[4];
   const double f = a_[5];
   const double g = a_[6];
   const double h = a_[7];
   const double k = a_[8];
   
   const double A = (e*k-f*h);
   const double B = (f*g-d*k);
   const double C = (d*h-e*g);
   const double D = (c*h-b*k);
   const double E = (a*k-c*g);
   const double F = (b*g-a*h);
   const double G = (b*f-c*e);
   const double H = (c*d-a*f);
   const double K = (a*e-b*d);
   const double Z = a*(e*k-f*h) + b*(f*g-k*d) + c*(d*h-e*g);

   if (Z != 0 ) {
      a_[0] = A;
      a_[1] = D;
      a_[2] = G;
      a_[3] = B;
      a_[4] = E;
      a_[5] = H;
      a_[6] = C;
      a_[7] = F;
      a_[8] = K;
      a_ = a_/Z;
   } else {
      a_[0] = 1.0;
      a_[1] = 0.0;
      a_[2] = 0.0;
      a_[3] = 0.0;
      a_[4] = 1.0;
      a_[5] = 0.0;
      a_[6] = 0.0;
      a_[7] = 0.0;
      a_[8] = 1.0;
   }
#endif

   /*
   gsl_matrix *a_ = gsl_matrix_alloc(NDIM,NDIM);
   gsl_matrix *b_ = gsl_matrix_alloc(NDIM,NDIM);
   for (int i=0;i<NDIM;i++) {
      for (int j=0;j<NDIM;j++) {
         gsl_matrix_set(a_,i,j,a[i*NDIM+j]);
      }
   }

   gsl_permutation *p = gsl_permutation_alloc(NDIM);
   int signum;

   gsl_linalg_LU_decomp (a_, p, &signum);
   gsl_linalg_LU_invert (a_, p, b_);

   //gsl_linalg_cholesky_decomp(a_);
   //gsl_linalg_cholesky_invert(a_);

   for (int i=0;i<NDIM;i++) {
      for (int j=0;j<NDIM;j++) {
         a[i*NDIM+j] = gsl_matrix_get(b_,i,j);
      }
   }
   gsl_matrix_free(a_);
   gsl_matrix_free(b_);
   gsl_permutation_free(p);
   */
}


double len2(vect inV); 
double len(vect inV); 
//double len2(vectInt inV); 
//double len(vectInt inV); 


#endif
