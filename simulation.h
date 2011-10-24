#ifndef SIMULATION_H
#define SIMULATION_H

#include "sph.h"
#include "vect.h"
#include "data.h"
#include "particle.h"

#define NBODIES 0

class Cbody {
   public:
      vect r;
      double mass;
};

class Csimulation {
	public:
      double maxdt;

      Csimulation(Cdata *_data); 
		void run();

      static Cbody bodies[NBODIES];
      static void calc_gravity_from_point_masses(Cparticle &p);
	private:
		void before_start();
		void before_middle();
		void before_end();

      Cdata *data;
      Csph *sph;
      
};

#endif
