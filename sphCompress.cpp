#include "sphCompress.h"

using namespace Nsph;

int CsphCompress::num_recalc_dens;
double CsphCompress::dt = 0;
double CsphCompress::hdt=0;
double CsphCompress::maxdt=0;
double CsphCompress::dtsig;

void CsphCompress::start() {
   data->newTimeStep();
	//start of verlet symplectic
	data->traverse<drift>();
}

void CsphCompress::middle() {
   data->traverse<init_hErr>();
	data->traverse<calc_h>();
	for (int i=0;i<NUM_POINT_ITERATIONS;i++) {
      data->traverse<init_density,if_hErr_greater_than>();
		data->neighbours<calc_density,if_hErr_greater_than>();
      num_recalc_dens = 0;
      data->traverse<update_hErr_and_calc_h,if_hErr_greater_than>();
      printf("Iteration %d, Number of particles to recalculate = %d\n",i,num_recalc_dens);
      if (num_recalc_dens == 0) break;
	}
   data->traverse<init_aom>();
   data->neighbours<calc_aom>();
	data->traverse<calc_pressure_spsound_and_pdr2>();	
   data->traverse<init_press_visc_force>();

   dtsig = 1000;
   
	data->neighbours<calc_press_visc_force>();

   dt = min(0.5*dtsig,maxdt);
   hdt = 0.5*dt;

	//finish verlet symplectic
	data->traverse<kick>();

   data->traverse<calc_dhdt_and_dalphdt>();
}

void CsphCompress::end() {
   data->traverse<drift_no_v>();
}

void CsphCompress::drift(Cparticle &p) {

	p.r += hdt*p.v;
   p.v0 = p.v;
	p.v += hdt*p.f;
   p.dens += hdt*p.dddt;
   p.u += hdt*p.dudt;
   p.h += hdt*p.dhdt;
   p.alpha += hdt*p.dalphdt;
}

void CsphCompress::drift_no_v(Cparticle &p) {

	p.r += hdt*p.v;
   p.dens += hdt*p.dddt;
   p.u += hdt*p.dudt;
   p.h += hdt*p.dhdt;
   p.alpha += hdt*p.dalphdt;
}
	
void CsphCompress::kick(Cparticle &p) {
   p.f = p.fp+p.fv+p.fg;
	p.v = p.v0 + dt*p.f;
}

void CsphCompress::calc_pressure_spsound_and_pdr2(Cparticle &p) {
   p.press = p.u*(GAMMA-1)*p.dens;
   p.spsound = sqrt(p.u*(GAMMA-1));
   p.pdr2 = p.press/(p.dens*p.dens);
}
   
bool CsphCompress::if_hErr_greater_than(Cparticle &p) {
   return (p.hErr > MAX_H_ERROR);
}
   
void CsphCompress::init_hErr(Cparticle &p) {
   p.hErr = 1;
}


void CsphCompress::update_hErr_and_calc_h(Cparticle &p) {
   double old_h = p.h;
   calc_h(p);
   p.hErr = abs((p.h-old_h)/old_h);
   if (p.hErr > MAX_H_ERROR) {
      num_recalc_dens++;
   }
}

void CsphCompress::init_press_visc_force(Cparticle &p) {
   p.fp = 0.0;
   p.fv = 0.0;

   p.dddt = 0.0;
   p.dudt = 0.0;

   p.maxvsig = 0.0;
}


void CsphCompress::calc_press_visc_force(Cparticle &pa, Cparticle &pb) {

	double r2 = len2(pa.r-pb.r);
	double r = sqrt(r2);
	double qa = r/pa.h;
	double qb = r/pb.h;
	vect dx = pa.r-pb.r;
	vect dv = pa.v-pb.v;
	double vdr = dot(dx,dv);
   double dr = 0.0;
   if (r!=0.0) dr = 1/r;
	double viss = vdr*dr;
   //cout << "abs(viss) = "<<abs(viss)<<" avspsound = "<<0.5*(pa.spsound+pb.spsound)<<endl;
	double vsig = pa.spsound + pb.spsound + 2*abs(viss);
   if (vsig > pa.maxvsig) pa.maxvsig = vsig;
   double visc = -viss*vsig*0.5*(pa.alpha+pb.alpha)/(pa.dens+pb.dens);
   dtsig = min(dtsig,min(pa.h,pb.h)/vsig);
	double Fa = F(qa,pa.h);
	double Fb = F(qb,pb.h);
	double dwp = pa.pdr2*Fa/pa.aom + pb.pdr2*Fb/pb.aom;
	double dwv = 0.5*visc*(Fa+Fb);
	pa.fp -= dx*pb.mass*dwp;
	pa.fv -= dx*pb.mass*dwv;

   double dddtInc = pb.mass*vdr*Fa/pa.aom;
   pa.dddt += dddtInc;
   pa.dudt += pa.pdr2*dddtInc;
}

