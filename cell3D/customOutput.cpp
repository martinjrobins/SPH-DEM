#include "customOutput.h"
#include "dataLL.impl.h"
#include "sph.h"
#include "misc.h"

using namespace Nsph;

inline void findCOMs(Cparticle &p, CglobalVars &g, vect* &com) {
   for (int i=0;i<5;i++) {
      if (len2(p.vhat)<pow(VREF/pow(10,i+1),2)) {
         const double r = sqrt(p.r[0]*p.r[0] + p.r[1]*p.r[1]);
         com[i+5][0] += r;
         com[i+5][1] += p.r[2];
         com[i+5][2] += 1.0;
      }
      if (p.shepSum<(0.3+i/10.0)) {
         const double r = sqrt(p.r[0]*p.r[0] + p.r[1]*p.r[1]);
         com[i][0] += r;
         com[i][1] += p.r[2];
         com[i][2] += 1.0;
      }
   }
}


void CcustomOutput::calcOutput(int outstep,CcustomSim *custSim,Cio_data_vtk *io) {
   fo << data->globals.time;
   //vector<vect> com(10);
   vect *com = (vect *)malloc(sizeof(vect)*10);;
   for (int i=0;i<10;i++) {
      com[i] = 0.0;
   }
   data->traverse<vect*&,findCOMs,ifDem>(com);
   for (int i=0;i<10;i++) {
      com[i][0] /= com[i][2];
      com[i][1] /= com[i][2];
      fo <<' '<<com[i][0]<<' '<<com[i][1];
   }
   fo << endl;
}





     
