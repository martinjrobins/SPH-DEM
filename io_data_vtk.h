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


/*
 * upper level I/O class. Contains globals and domain I/O functionality, as
 * well as adding the ability to read/write particle data using the VTK file format
 */
class Cio_data_vtk: public CioGlobals, public CioDomain {
public:
      /*
       * constructor sets filename
       */
      Cio_data_vtk(string _filename,CglobalVars *g): CioGlobals(_filename,g),CioDomain(_filename,g) {
         filename = _filename;
      };

      /*
       * read/write restart datafile (used to restart a simulation from scratch)
       */
      void readRestart(int timestep,particleContainer *outPs,CglobalVars *globals);
      void writeRestart(int timestep,particleContainer &ps,CglobalVars *globals);

      /*
       * read/write regular output file
       */
      void readOutput(int timestep,particleContainer *outPs,CglobalVars *globals);
      void writeOutput(int timestep,particleContainer &ps,CcustomSimBase &customSim,CglobalVars *globals);

      /*
       * read output file using CSIRO format
       */
      void readCSIRO(int timestep,particleContainer *outPs,CglobalVars *globals);

      void writeAux(int timestep,particleContainer &ps,vector<double> theData,const char *name,CglobalVars *globals);
      void writeAuxNoData(int timestep,particleContainer &ps,const char *name,CglobalVars *globals);
      void writeGrid(int timestep,particleContainer &ps,vectInt &gridDims,CglobalVars *globals);

      /*
       * set/get base filename. 
       * filenames for regular output datafile is:
       * "baseNameNNNNNNN.pvtu" where NNNNNNN is the output timestep number
       * and baseName is the base filename
       *
       * for a restart file, the total filename is:
       * "baseNameRestartNNNNNNN.pvtu"
       */
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
