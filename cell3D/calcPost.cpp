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
}

inline void calcPorosityAndVelocity(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g) {
	p.v = 0.0;
	p.dens = 0.0;
	p.porosity = 0.0;
	p.shepSum = 0.0;
        double sum = 0.0;
	int n = neighbrs.size();
        //if (n>0) cout<<n<<endl; 
	for (int i=0;i<n;i++) {
	   Cparticle *pn = neighbrs[i];
	   if (pn->iam==boundary) continue;
	   const vect dx = p.r-pn->r;
	   const double r = len(dx);
	   const double dvol = pn->mass/pn->dens;
	   const double hav = pn->h;
	   const double q = r/hav;
	   const double Wab = W(q,hav);
           sum += dvol*Wab;
	   p.v += pn->v*dvol*Wab;
	   p.dens += pn->mass*Wab;
	   p.porosity += pn->porosity*dvol*Wab;
	}
	p.shepSum = sum;
	if (sum != 0.0) {
	   p.v /= sum;
	   p.dens /= sum;
	   p.porosity /= sum;
	}
	//cout <<p.r<<" "<<p.h<<" "<<p.porosity<<endl;
}

inline void findCOMs(Cparticle &p, CglobalVars &g, vect* &com) {
   for (int i=0;i<5;i++) {
      if (len2(p.v)<pow(VREF/pow(10,i+1),2)) {
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



int main(int argc, char *argv[]) {
   if (argc != 4) {
      cout << "Usage: averageVelocityField infilename startStep endStep" << endl;
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


   ofstream fo;
   fo.open("com.dat");
   fo << "time dry0.3:r,h dry0.4 dry0.5 dry0.6 dry0.7 still10 still100 still1000 still10000 still100000"<<endl;

   Nmisc::setupEntireDomain(globals.procDomain,globals.procNeighbrs);


   cout << "Processor "<<globals.mpiRank<<" of "<<globals.mpiSize<<" has domain ";
   for (int i=0;i<NDIM*2;i++) {
      cout<<globals.procDomain[i]<<' ';
   }
   cout << " and neighbrs "<<globals.procNeighbrs<<endl;
  
   cout <<"read first data set"<<endl;
   particleContainer ps;
   io_data.readOutput(startStep,&ps,&globals);
   cout <<"creating data structure"<<endl;
   CdataLL *data= new CdataLL(ps,globals,false);


   particleContainer gridPs;
   vectInt gridDims;
   data->functionOverGrid<init>(gridPs,BMIN,BMAX,gridDims);
   //CdataLL *data_grid= new CdataLL(gridPs,globals,false);

   //io_data.setFilename(outName,&globals);
   for (int nstep=startStep;nstep<=endStep;nstep++) {
      cout << "Timestep "<<nstep<<" Time "<<globals.time<<" Post-Process "<<double(nstep-startStep)/double(endStep-startStep)*100<<"\% complete"<<endl;
      globals.outstep = nstep;
      //globals.time = nstep*DT;
      cout <<"reading globals"<<endl;
      io_data.readGlobals(nstep,&globals);
      cout <<"reading output"<<endl;
      io_data.readOutput(nstep,&ps,&globals);
      //data->traverse<setMass>();
      cout <<"updating particle"<<endl;
      data->updateParticles();
      data->reset();
      data->neighboursUsing<calcPorosityAndVelocity>(gridPs);

      //io_data.writeOutput(globals.outstep,gridPs,&globals);
      io_data.writeGrid(globals.outstep,gridPs,gridDims,&globals);

      fo << globals.time;
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
   fo.close();
}

