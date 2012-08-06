#include "customConstants.h"
#include "vect.h"
#include "sphIncompress.h"
#include "dataLL.impl.h"
#include "io_data_vtk.h"
#include "particle.h"
#include "misc.h"
#include <cstdlib>
#include <string>

struct test {
   public:
      test() {};
static double readEval1(Cparticle &p) {
   return p.eval1;
}
static void writeEval1(Cparticle &p,double eval1) {
   p.eval1 = eval1;
}       
};

inline void calcGradEval1(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g,CdataLL &d2) {
   //calc gradient of largest eigenvalue
   p.grad_eval1 = 0;
   const int n = neighbrs.size();
   for (int i;i<n;i++) {
      Cparticle *pn = neighbrs[i];
      const vect dr = p.r-pn->r;
      const vect dEval1 = p.eval1-pn->eval1;
      const double r = len(dr);

#ifdef VAR_H_CORRECTION
      const vect gradWa = gradW(p,*pn,dr,r/p.h,p.h);
#else
      const vect gradWa = 0.5*(gradW(p,*pn,dr,r/p.h,p.h)+gradW(p,*pn,dr,r/pn->h,pn->h));
#endif
      p.grad_eval1 += pn->mass*dot(dEval1,gradWa);
   }
   p.grad_eval1 *= 1.0/p.dens;
}



inline void calcHallerLCSExtractParameters(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g,CdataLL &d2) {
   calcGradEval1(p,neighbrs,g,d2);
   p.Z_surface = dot(p.grad_eval1,p.evec1);
   cout <<p.Z_surface<<" "<<p.grad_eval1<<" "<<p.evec1<<endl;
}



int main(int argc, char *argv[]) {
   if (argc != 4) {
      cout << "Usage: post infilename startStep endStep" << endl;
      return(-1);
   }
   string filename = argv[1];
   int startStep = atoi(argv[2]);
   int endStep = atoi(argv[3]);
   CglobalVars globals;
   Cio_data_vtk io_data(filename,&globals);

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &(globals.mpiSize));
   MPI_Comm_rank(MPI_COMM_WORLD, &(globals.mpiRank));

   Nmisc::setupEntireDomain(globals.procDomain,globals.procNeighbrs);

   cout << "Processor "<<globals.mpiRank<<" of "<<globals.mpiSize<<" has domain ";
   for (int i=0;i<NDIM*2;i++) {
      cout<<globals.procDomain[i]<<' ';
   }
   cout << " and neighbrs "<<globals.procNeighbrs<<endl;
  
   cout <<"read first data set"<<endl;
   particleContainer ps_start;
   io_data.readCSIRO(startStep,&ps_start,&globals);
   cout <<"creating data structure"<<endl;
   CdataLL *data_start = new CdataLL(ps_start,globals,false);


   CglobalVars globals2;
   Cio_data_vtk io_data2(filename,&globals2);
   particleContainer ps_end;
   io_data2.readCSIRO(endStep,&ps_end,&globals2);
   CdataLL *data_end = new CdataLL(ps_end,globals2,false);
   
   globals.outstep = 1;

   cout <<"calculating Right Cauchy-Green Deformation Tensor..."<<endl;
   data_start->neighboursGroup<CdataLL &,Nmisc::calcRightCauchyGreenTensor,Nsph::ifSph>(*data_end);
   cout <<"syncing...."<<endl;
   data_start->syncParticlesBetweenProcs<double,test::readEval1,test::writeEval1>();
   cout <<"calculating Haller LCS Parameters..."<<endl;
   data_start->neighboursGroup<CdataLL &,calcHallerLCSExtractParameters,ifSph>(*data_end);

   string outName;
   outName = "HallerLCS_ForwardTime";

   io_data.setFilename(outName,&globals);
   io_data.writeOutput(globals.outstep,ps_start,&globals);
}

