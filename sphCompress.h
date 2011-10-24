#ifndef SPHCOMPRESS_H
#define SPHCOMPRESS_H

#include "data.h"
#include "vect.h"
#include "particle.h"
#include "sph.h"

#define MAX_H_ERROR 0.01

#define NUM_POINT_ITERATIONS 6

class CsphCompress {
   public:
      static int num_recalc_dens;
      static double dt,hdt,maxdt;
      static double dtsig;

      CsphCompress(Cdata *_data) {
         data = _data;
      };
      
      void start();
      void middle();
      void end();

      static void kick(Cparticle &p);
      static void drift(Cparticle &p);
      static void drift_no_v(Cparticle &p);

      static void calc_pressure_spsound_and_pdr2(Cparticle &p);
      static void init_hErr(Cparticle &p);
      
      static bool if_hErr_greater_than(Cparticle &p);
      static void update_hErr_and_calc_h(Cparticle &p);

      static void init_press_visc_force(Cparticle &p);
      static void calc_press_visc_force(Cparticle &pa,Cparticle &pb);
      

   private:
      Cdata *data;

};

#endif
