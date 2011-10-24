void simulation::run() {
	//go throught timesteps
	shouldstop = time > MAXTIME || nstep > MAXSTEP;
	while (!shouldstop) {
		//do stuff before timestep starts
		before_start();
		//do the start of the timestep
		gsph.start();
		//do stuff before the middle of the timestep
		before_middle();
		//do the middle of the timestep
		gsph.middle();
		//do stuffbefore the end of the timestep
		before_end();
		//do the end of the timestep
		gsph.end();
	}
}

void simulation::before_start() {
   if (nstep%IOSTEP==0) {
	   io_data.write();
   }
   nstep++;
}

void simulation::before_middle() {
   time += gsph.hdt;
}

void simulation::before_end() {
   time += gsph.hdt;
}
