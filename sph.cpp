#include "sph.h"


void Nsph::init_aom(Cparticle &p) {
   p.aom= 1;
   calc_aom(p,p);
}


void Nsph::calc_aom(Cparticle &pa,Cparticle &pb) {

   double qa = len(pa.r-pb.r)/pa.h;
   pa.aom -= (1.0/(NDIM*pow(HFAC,NDIM)*pa.mass))*pb.mass*(dKdq(qa,pa.h)*qa+NDIM*K(qa,pa.h));
   //cout << "aom increment = "<<(1.0/(NDIM*pow(HFAC,NDIM)*pa.mass))*pb.mass*(Csph::dKdq(qa,pa.h)*qa+NDIM*Csph::K(qa,pa.h))<<endl;
   //cout << "right most bit = "<<(Csph::dKdq(qa,pa.h)*qa+NDIM*Csph::K(qa,pa.h))<<endl;
   //pa.aom -= (1/NDIM)*(1/pow(HFAC,NDIM))*(pb.mass/pa.mass)*(Csph::dKdq(qa,pa.h)*qa+NDIM*Csph::K(qa,pa.h));
   //pa.aom = 1;
}
   



void Nsph::calc_h(Cparticle &p) {
   p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
   if (p.h < 0) cout << "h is crap now!!!! p.dens= "<<p.dens<<endl;
}

	
void Nsph::zero_density(Cparticle &p) {
   p.dens = 0.0;
}

void Nsph::init_density(Cparticle &p) {
   p.dens = 0.0;
   calc_density(p,p);
}

void Nsph::calc_density(Cparticle &pa, Cparticle &pb) {
   
	double r = len(pa.r-pb.r);
	pa.dens += pb.mass*W(r/pa.h,pa.h);
}


void Nsph::calc_dhdt_and_dalphdt(Cparticle &p) {
   double graddotv = -p.dddt/p.dens;
   p.dhdt = 0.5*p.h*graddotv;
   double tou = p.h/(0.1*p.maxvsig);
   double s = max(-graddotv,0.0)*(2.0-p.alpha);
   p.dalphdt = -(p.alpha-MIN_ALPHA)/tou + s;
}

double Nsph::F(double q, double h) {
   if (q<=1.0) {
      return (1/pow(h,NDIM+2))*WCON*(-2.0+ 1.5*q);
   }
   else if (q<=2.0) {
      return -(1/pow(h,NDIM+2))*3.0*(WCON/6.0)*pow(2.0-q,2)/q;
   }
   else {    
      return 0.0; 
   }
}

double Nsph::dKdq(double q, double h) {
   if (q<=1.0) {
      double q2 = pow(q,2);
      return WCON*(-2.0*q + 1.5*q2);
   }
   else if (q<=2.0) {
      return -3.0*(WCON/6.0)*pow(2.0-q,2);
   }
   else {    
      return 0.0; 
   }
}

double Nsph::W(double q, double h) {
   return (1/pow(h,NDIM))*K(q,h);
}

double Nsph::K(double q, double h) {
   if (q<=1.0) {
      double q2 = pow(q,2);
      double q3 = q*q2;
      return WCON*(2.0/3.0 - q2 + 0.5*q3);
   }
   else if (q<=2.0) {
      return (WCON/6.0)*pow(2.0-q,3);
   }
   else {    
      return 0.0;
   }
}
