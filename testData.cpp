
#include <iostream>
#include <algorithm>
#include <vector>
#include "vect.h"
#include "particle.h"
#include "data.h"
using namespace std;

void traverseFunct(Cparticle &ppa) {
   cout << "Traversing through particle with tag = " << ppa.tag << endl;
}

void neighbrsFunct(Cparticle &ppa, Cparticle &ppb) {
   cout << "particle " << ppa.tag << " and particle " << ppb.tag << " are neighbours" << endl;
}
int main() {
   cout << "Creating data object..." << endl;
   Cdata *thedata;
   vector<Cparticle> inParticles;
   inParticles.push_back(Cparticle(vect(-100,-123)));
   inParticles.push_back(Cparticle(vect(1,0.2)));
   inParticles.push_back(Cparticle(vect(1,0.1)));
   inParticles.push_back(Cparticle(vect(1,1)));
   
   for (int i=0;i<4;i++) {
      inParticles[i].mass = 0.25;
      inParticles[i].h = 0.3;
      inParticles[i].tag = i;
   }
   
   cout << "Constructing data structure..." << endl;
   thedata = new Cdata(inParticles);

   cout << "Calling traverse function..." << endl;
   thedata->traverse<traverseFunct>();

   cout << "Calling neighbours function..." << endl;
   thedata->neighbours<neighbrsFunct>();

   cout << "Calling neighbours function again..." << endl;
   thedata->neighbours<neighbrsFunct>();

   cout << "reconstruct" << endl;
   thedata->newTimeStep();
   
   cout << "Calling traverse function..." << endl;
   thedata->traverse<traverseFunct>();

   cout << "Calling neighbours function..." << endl;
   thedata->neighbours<neighbrsFunct>();

   cout << "Calling neighbours function again..." << endl;
   thedata->neighbours<neighbrsFunct>();


   cout << "All done. Have a nice day" << endl;
   return 0;
}
   
   

   
