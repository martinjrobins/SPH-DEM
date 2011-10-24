Csph::start() {
	//start of verlet symplectic
	gdata.traverse<drift>();
}

Csph::middle() {
	for (int i=0;i<NO_POINT_ITERATIONS;i++) {
      gdata.neighbours<zero_density>();
		gdata.neighbours<calc_density>();
	   sumherror2 = 0;	
		no_big_herror = 0;
		gdata.traverse<calc_h>();
		//TODO:output point iteration errors here
	}
	gdata.traverse<calc_pressure_and_aom>();	
   gdata.traverse<zero_forces>();
	gdata.traverse<calc_gravity_from_point_masses>();
	gdata.neighbours<calc_press_visc_force>();
	gdata.traverse<calc_timestep>();

	//finish verlet symplectic
	gdata.traverse<kick>();
}

Csph::end() {
}

void Csph::drift(particle *pp)
	particle->x += hdt*particle->v;
}

	
void Csph::kick(particle *pp) {
	pp->v += dt*pp->force;
	pp->v *= DAMP;
	pp->x += hdt*pp->v;
}

void Csph::calc_h(particle *pp) {
	real hold = pp->h;
	pp->h = HFAC*sqrt(pp->mass/pp->dens);
	real herror = (hold-pp->h)**2/PSEP**2;
	sumherror2 += herror;
	if (herror>BIG_H_ERROR) no_big_herror++;
}

void Csph::zero_forces(particle *pp) {
   pp->fg = 0.0;
   pp->fp = 0.0;
   pp->fv = 0.0;
}

void Csph::calc_gravity_from_point_masses(particle *pp) {
   for (int i=0;i<gsimulation->nbodies;i++) {
      vect bx = gsimulation->bodies[i]->x;
      real bm = gsimulation->bodies[i]->mass;
      real r2 = len2(pp->x-gsimulation->bodies[i]->x);
      real r1 = sqrt(r2);
      real r3 = r1*r3;
      vect dx = bx-pp->x;
   
      pp->fg +=  dx*bm/(r3+GFAC);
   }
}
	
void Csph::zero_density(particle *pp) {
   pp->dens = 0.0;
}

void Csph::density(particle *ppa, particle *ppb) {
	real r2 = len2(ppa->x-ppb->x);
	ppa->dens += ppa->mass*W(r,r2,ppa->h);
}

void Csph::calc_press_visc_force(particle *ppa, particle *ppb) {
	real r2 = len2(ppa->x-ppb->x);
	real r = sqrt(r2);
	real qa = r/ppa->h;
	real qb = r/ppb->h;
	vect dx = ppa->x-ppb->x;
	vect dv = ppa->v-ppb->v;
	real vdr = dx.dot(dv);
	real viss = vdr/r;
	real vsig = ppa->spsound + ppb->spsound + 2*abs(viss);
	dtsig = min(dtsig,min(ppa->h,ppb->h)/vsig);
	real Fa = F(qa,ha);
	real Fb = F(qb,hb);
	real dwp = ppa->pdr2*Fa/ppa->aom + ppb->pdr2*Fb/ppb->aom;
	real dwv = 0.5*visc*(Fa+Fb);
	ppa->fp += dx*ppa->mass*dwp;
	ppa->fv += dx*ppa->mass*dwv;
}


real Csph::F(real q, real h) {
   if (q<=1.0) {
      return (1/h)^(NDIM+2)*WCONA*(-2.0+ 1.5*q);
   }
   else if (q<=2.0) {
      return -(1/h)^(NDIM+2)*3.0*WCONB*(2.0-q)^2/q;
   }
   else {    
      return 0.0; 
   }
}

real Csph::W(real q, real h) {
   if (q<=1.0) {
      real q2 = q**2;
      real q3 = q*q2;
      return (1/h)^(NDIM+2)*WCONA*(2.0/3.0 - q2 + 0.5*q3);
   }
   else if (q<=2.0) {
      return (1/h)^(NDIM+2)*WCONB*(2.0-q)^3;
   }
   else {    
      return 0.0 
   }
}

