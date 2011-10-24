#ifndef IO_DATA_VTK_H
#define IO_DATA_VTK_H

#include "customConstants.h"
#include "particle.h"
#include "globalVars.h"
#include "customSimBase.h"
#include "ioGlobals.h"
#include "ioDomain.h"
#include <vector>
#include <string>


class Cio_data_vtk: public CioGlobals, public CioDomain {
	public:
      Cio_data_vtk(string _filename,CglobalVars *g): CioGlobals(_filename,g),CioDomain(_filename,g) {
         filename = _filename;
      };
      void readRestart(int timestep,particleContainer *outPs,CglobalVars *globals);
      void readOutput(int timestep,particleContainer *outPs,CglobalVars *globals);
      void readCSIRO(int timestep,particleContainer *outPs,CglobalVars *globals);
      void writeOutput(int timestep,particleContainer &ps,CcustomSimBase &customSim,CglobalVars *globals);
      void writeAux(int timestep,particleContainer &ps,vector<double> theData,const char *name,CglobalVars *globals);
      void writeAuxNoData(int timestep,particleContainer &ps,const char *name,CglobalVars *globals);
      void writeRestart(int timestep,particleContainer &ps,CglobalVars *globals);
      void writeGrid(int timestep,particleContainer &ps,vectInt &gridDims,CglobalVars *globals);
		void setFilename(string _filename,CglobalVars *g) { 
         CioGlobals::setFilename(_filename,g);
         CioDomain::setFilename(_filename,g);
         filename = _filename; 
      };
      string getFilename() { return filename; }
   private:
      string filename;
      char *dimNames[3]; 
      char *vdimNames[3];
};


#endif
