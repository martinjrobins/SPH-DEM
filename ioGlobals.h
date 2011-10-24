#ifndef IOGLOBALS_H
#define IOGLOBALS_H

#include "particle.h"
#include <string>
#include <mpi.h>
#include <fstream>
#include "globalVars.h"

/*
 * I/O class for global variables. 
 * The global variables are stored in a space separated column list in an ascii
 * file.
 * each newline of the file is a new output timestep.
 */
class CioGlobals {
public:
   /*
    * set the base filename. The output filename will be 
    * "baseNameGlobalsN.dat" where "baseName" is the base filename chosen and
    * "N" is the cpu number. cpu 0 has the overall global variables (for all
    * the cpus)
    */
   void setFilename(string _filename,CglobalVars *g);

   /*
    * constructor sets base filename
    */
   CioGlobals(string _filename,CglobalVars *g) {
      setFilename(_filename,g);
   };

   /*
    * read/write to a given output timestep
    */
   void readGlobals(int outStep,CglobalVars *g);
   void writeGlobals(int outStep,CglobalVars *g);

   /*
    * scan globals file for a given time
    */
   int findClosestStepToTime(double time);
private:
   string filename;
   fstream fio;
   int myGetLine(CglobalVars *g);
   void myWriteLine(int currStep,CglobalVars *g);
   void openAndSeekToStep(int outStep,CglobalVars *g);
};

#endif
