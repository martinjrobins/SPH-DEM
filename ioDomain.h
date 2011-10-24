#ifndef IODOMAIN_H
#define IODOMAIN_H

#include <string>
#include <fstream>
#include "globalVars.h"

/* 
 * I/O class for the domain information (extent of the simulated domain for
 * the given CPU)
 */
class CioDomain {
   public:
      void setFilename(string _filename,CglobalVars *g);

      /*
       * constructor sets the base filename
       * the total filename for the domain info is:
       * "baseNameDomainN.dat" where N is the CPU number
       */
      CioDomain(string _filename,CglobalVars *g) {
         setFilename(_filename,g);
      }

      /*
       * read/write domain info for a given output timestep. Note that at the
       * moment outStep does nothing, but in the future the domain will be
       * allowed to change over time
       */
      void readDomain(int outStep,CglobalVars *g);
      void writeDomain(int outStep,CglobalVars *g);

      /*
       * depreciated functions
       */
      void readDomainOld(int outStep,CglobalVars *g);
      void writeDomainOld(int outStep,vector<vector<double> > vprocDomain, vector<Array<int,NDIM> > vprocNeighbrs);

   private:
      string filename;
      fstream fio;
};

#endif
