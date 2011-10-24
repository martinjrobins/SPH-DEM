#include "simulation.h"


Cbody Csimulation::bodies[NBODIES];
         

Csimulation(Cdata _data) {
         maxdt=1000;
         data = _data;
         sph = new Csph(data);
         sph.maxdt = MAXTIME/MAXSTEP;
};

void simulation::run() {
	//go throught timesteps
	shouldstop = time > MAXTIME || nstep > MAXSTEP;
	while (!shouldstop) {
		//do stuff before timestep starts
		before_start();
		//do the start of the timestep
		sph.start();
		//do stuff before the middle of the timestep
		before_middle();
		//do the middle of the timestep
		sph.middle();
		//do stuffbefore the end of the timestep
		before_end();
		//do the end of the timestep
		sph.end();
	}
}

void simulation::before_start() {
   if (nstep%IOSTEP==0) {
	   io_data.write(nstep);
   }
   nstep++;
}

void simulation::before_middle() {
   time += gsph.hdt;
   gsph->data->traverse<calc_gravity_from_point_masses>();
}

void simulation::before_end() {
   time += gsph.hdt;
}


static void calc_gravity_from_point_masses(Cparticle &p) {

   for (int i=0;i<NBODIES;i++) {
      vect bx = bodies[i].r;
      double bm = bodies[i].mass;
      double r2 = len2(p.x-bodies[i].r);
      double r1 = sqrt(r2);
      double r3 = r1*r3;
      vect dx = bx-p.r;
   
      p.fg +=  dx*bm/(r3+GFAC);
   }
}
