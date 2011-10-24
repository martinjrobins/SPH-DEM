#include "constants.h"
#include "particle.h"
#include "io_data.h"
#include <vector>
#include <iostream>

int main() {
   vector<Cparticle> ps;
   Cparticle p;
   for (int i=0;i<NX;i++) {
      for (int j=0;j<=NY;j++) {
         p.r = (i+0.5)*PSEP+RMIN[0],j*PSEP+RMIN[1];
         p.dens = DENS;
         p.mass = PSEP*PSEP*DENS;
         p.h = H;
         if (j==0) {
            p.v = WALLSPEED,0;
            p.iam = boundary;
            p.norm = 0.0,1.0;
         } else if (j==NY) {
            p.v = -WALLSPEED,0;
            p.iam = boundary;
            p.norm = 0.0,-1.0;
         } else {
            p.v = 0,0;
            p.iam = sph;
         }
         p.alpha = ALPHA;
         ps.push_back(p);
      }
   }

   cout << "Total number of particles = " << ps.size() << endl;

   Cio_data ioFile("couette");
	ioFile.write(1,ps);

}

