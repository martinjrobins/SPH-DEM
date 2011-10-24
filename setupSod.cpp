#include "constants.h"
#include "particle.h"
#include "io_data.h"
#include <vector>
#include <iostream>

int main() {
	
   const double psep = 1.0/NX;
   const double dr = 0.125; //density on the right
   const double dl = 1.0;   //  density on the far left
   const double gam1 = GAMMA-1.0;
   const double psepr = psep;
   const double psepl = psep*dr/dl;
   const double xmin = -2.0;
   const double xmax = 2.0;
   const double xminstill = xmin+4*HFAC*psepr;
   const double xmaxstill = xmax-4*HFAC*psepr;
   const double ur = 0.1/(dr*gam1);    //  utherm  on  right
   const double ul = 1.0/(dl*gam1);    // utherm  on  left
   const double al = 0.5*psepr;
   const double pmref = dr*psepr;

   cout << "adding first particle" << endl;
   vector<Cparticle> ps;
   Cparticle p;
   cout << "adding first particle" << endl;
   p.r[0] = xmin;
   p.dens = dl;
   p.mass = dr*psepr;
   p.h = HFAC*psepl;
   p.u = ul;
   p.v[0] = 0.0;
   p.alpha = 0.5;
   p.iam = immovable;
   p.colour = 2;
   p.tag = 1;
   cout << "adding first particle" << endl;
   ps.push_back(p);
   cout << "adding first particle" << endl;
   cout << ps.size();
   cout << "adding second particle" << endl;
  
   cout << "adding second particle" << endl;
   p.r[0] = xmin+psepl;
   p.dens = dl;
   p.mass = dr*psepr;
   p.h = HFAC*psepl;
   p.u = ul;
   p.iam = immovable;
   p.colour = 2;
   p.v[0] = 0.0;
   p.alpha = 0.5;
   p.tag = 2;
   cout << "adding second particle" << endl;
   ps.push_back(p);

   cout << "adding the rest" << endl;
   while(ps.back().r[0] < xmax) {
      cout << "there is currently " << ps.size() << "particles in the ps" << endl;
      p.r[0] = ps[ps.size()-2].r[0] + 2.0*pmref/ps[ps.size()-1].dens;
      cout << "added particle at x = " << p.r[0] << endl;
      if ((p.r[0]<xminstill) || (p.r[0]>xmaxstill)) {
        p.iam = immovable;
        p.colour = 2;
      } else {
        p.iam = sph;
        p.colour = 1;
      }
      double diff = p.r[0]/al;
      if (diff<-10) {
         p.dens = dl;
         p.u = ul;
         p.h = HFAC*psepl;
      } else if (diff>10) {
         p.dens = dr;
         p.u = ur;
         p.h = HFAC*psepr;
      } else {
         double exx = exp(diff);
         p.dens = (dl + dr*exx)/(1+exx);
         p.u  = (ul + ur*exx)/(1.+exx);
         p.h = HFAC*psep*dr/p.dens;
      }
      p.mass = dr*psepr;
      p.v[0] = 0.0;
      p.alpha = 0.5;
      p.tag = ps.size()+1;
      ps.push_back(p);
   }

   cout << "Total number of particles = " << ps.size() << endl;

   Cio_data ioFile("ssod.vtk");
   cout <<"boo"<<endl;
	ioFile.write(1,ps);

}

