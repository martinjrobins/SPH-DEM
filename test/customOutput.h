#ifndef CUSTOMOUTPUT_H
#define CUSTOMOUTPUT_H

//#include <fftw3.h>
#include <ext/slist>
#include "customOutputBase.h"
#include "customSim.h"
#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"
#include "misc.h"

namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later

using namespace Nmisc;

class CcustomOutput : public CcustomOutputBase {
   public:
      CcustomOutput(CdataLL *_data): CcustomOutputBase(_data) {
         //_data->globals.custom.resize(2);

         rng = gsl_rng_alloc(gsl_rng_ranlxd1);
         gsl_rng_set(rng,time(NULL));

         //fftw_in = (double *)fftw_malloc(sizeof(double)*(NX+1)*(NY+1));
         //fftw_out = (fftw_complex *)malloc(sizeof(fftw_complex)*(NX+1)*((int)floor((NY+1)/2)+1));
         //cout <<"Creating FFTW plan..."<<flush;
         //p = fftw_plan_dft_r2c_2d(NX+1,NY+1,fftw_in,fftw_out,FFTW_MEASURE|FFTW_DESTROY_INPUT);
         //cout <<"Finished."<<endl;
      }
      void calcOutput(int outstep,CcustomSim *custSim, Cio_data_vtk *io);
      void calcPostProcess(int outstep,CcustomSim *custSim, Cio_data_vtk *io);
   private:
      //pair doubling time data
      Sgi::slist<CpairData> pairData;
      vector<ofstream *> openFiles;
      gsl_rng *rng;

      //fft data
      //fftw_plan p;
      //fftw_complex *fftw_out;
      //double *fftw_in;

};


#endif
