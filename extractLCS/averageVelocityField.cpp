#include "customConstants.h"
#include "vect.h"
#include "sphIncompress.h"
#include "dataLL.impl.h"
#include "io_data_vtk.h"
#include "particle.h"
#include "misc.h"
#include <cstdlib>
#include <string>



inline void init(Cparticle &pa, Cparticle &pb,CglobalVars &g) {
   pa.tmp = 0;
}

inline void rotate(Cparticle &p, CglobalVars &g, double theta) {
   vect rotated;
   rotated[0] = cos(theta)*p.r[0] - sin(theta)*p.r[1];
   rotated[1] = sin(theta)*p.r[0] + cos(theta)*p.r[1];
   rotated[2] = p.r[2];
   p.r = rotated;
   rotated[0] = cos(theta)*p.v[0] - sin(theta)*p.v[1];
   rotated[1] = sin(theta)*p.v[0] + cos(theta)*p.v[1];
   rotated[2] = p.v[2];
   p.v = rotated;
   p.vhat = rotated;
}

inline void averageVelocity(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g) {
	vect vel = 0.0;
	int n = neighbrs.size();
	const double theta = g.time;
	for (int i=0;i<n;i++) {
	   Cparticle *pn = neighbrs[i];
	   if (pn->iam==boundary) continue;
	   const vect dx = p.r-pn->r;
	   const double r = len(dx);
	   const double dvol = pn->mass/pn->dens;
	   const double hav = pn->h;
	   const double q = r/hav;
	   const double Wab = W(q,hav);
	   vel += pn->vhat*dvol*Wab;
	}
	vect rotVel;
	rotVel[0] = -p.r[1];
	rotVel[1] = p.r[0];
	vel = vel - rotVel;
	p.vhat = p.tmp*p.vhat + vel;
	p.tmp = p.tmp + 1;
	p.vhat /= p.tmp;
	p.v = p.vhat;
}

int main(int argc, char *argv[]) {
   if (argc != 4) {
      cout << "Usage: averageVelocityField infilename startStep endStep dt" << endl;
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
   particleContainer ps;
   io_data.readCSIRO(startStep,&ps,&globals);
   cout <<"creating data structure"<<endl;
   CdataLL *data= new CdataLL(ps,globals,false);
   particleContainer gridPs;

   vectInt gridDims;
   data->functionOverGrid<init>(gridPs,BMIN,BMAX,gridDims);
   CdataLL *data_grid= new CdataLL(gridPs,globals,false);

   for (int nstep=startStep;nstep<=endStep;nstep++) {
      cout << "Timestep "<<nstep<<" Time "<<globals.time<<" Post-Process "<<double(nstep-startStep)/double(endStep-startStep)*100<<"\% complete"<<endl;
      globals.outstep = nstep;
      globals.time = nstep*DT;
      cout <<"reading globals"<<endl;
      io_data.readGlobals(nstep,&globals);
      cout <<"reading output"<<endl;
      io_data.readCSIRO(nstep,&ps,&globals);
      //data->traverse<setMass>();
      cout <<"updating particle"<<endl;
      data->updateParticles();
      data->reset();
      data_grid->traverse<rotate>(globals.time);
      data->neighboursUsing<averageVelocity>(gridPs);
      data_grid->traverse<rotate>(-globals.time);
   }
   string outName;
   outName = "averageVel";

   io_data.setFilename(outName,&globals);
   //io_data.writeOutput(globals.outstep,gridPs,&globals);
   io_data.writeGrid(globals.outstep,gridPs,gridDims,&globals);

}

