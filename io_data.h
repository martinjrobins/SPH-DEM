#include "particle.h"
#include "customConstants.h"
#include "customSim.h"
#include <vector>
#include <string>


class Cio_data {
	public:
      Cio_data(string _filename) {
         filename = _filename;
         dimNames[0] = "x";
         dimNames[1] = "y";
         dimNames[2] = "z";
         vdimNames[0] = "px";
         vdimNames[1] = "py";
         vdimNames[2] = "pz";
      };
		void read(int timestep,vector<Cparticle> *outPs,double *outTime,int *outNstep);
		void write(int timestep,vector<Cparticle> &ps,CcustomSim &customSim);
      void writeRestart(int timestep,vector<Cparticle> &ps);
      void writeGrid(int timestep,vector<Cparticle> &ps,vectInt &gridDims);
		void setFilename(string _filename) { filename = _filename; };
   private:
      string filename;
      char *dimNames[3]; 
      char *vdimNames[3];
};
