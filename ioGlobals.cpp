
#include "ioGlobals.h"

void CioGlobals::setFilename(string _filename,CglobalVars *g) {
   char str[20];
   sprintf(str,"Globals%d.dat",g->mpiRank);
   filename = _filename+str;
   fio.open(filename.c_str());

   if (!fio.is_open()) {
      cout << "Could not open Global Vars file: "<<filename<<". Create it..."<<endl;
      ofstream fo;
      fo.open(filename.c_str());
      if (!fo.is_open()) {
         cout <<"Crap! Still could not open Global Vars file, exiting..."<<endl;
         exit(-1);
      }
      fo <<"# OutStep SPHStep Time linMomX linMomY angMom KE eBForce eElast eElastExact eViscF eViscB eFF eTotal wtTotalStep wtTotalForceCalc wtTotalDddtCalc wtTotalReset wtTotalUpdateDomain wtTotalMPI maxV aveDens varDens aveDensFromMass n nSph maxF maxFF rmsFF edElastdt edViscFdt edFFdt"<<'\n';
      fo.close();
   }

   fio.close();
}

int CioGlobals::myGetLine(CglobalVars *g) {
   int currStep = -1;
   if (!fio.eof()) {
      string line;
      getline(fio,line);
      istringstream ss(line);

      //TODO: 2D specific!!!!
      double tmp;
      ss >> currStep >>g->sphStep >> g->time >> g->linMom[0]>>g->linMom[1]>>g->angMom>>g->eKE>>g->eBForce>>g->eElast>>g->eElastExact>>g->eViscF>>g->eViscB>>g->eFF>>g->eTotal>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>g->maxV>>g->aveDens>>g->varDens>>g->aveDensFromMass>>g->n>>g->nSph>>g->maxF>>g->maxFF>>g->rmsFF>>g->edElastdt>>g->edViscFdt>>g->edFFdt;
   }
   return currStep;
}

void CioGlobals::myWriteLine(int currStep,CglobalVars *g) {
   //TODO: 2D specific!!!!
   fio << currStep <<' '<< g->sphStep <<' '<<g->time<<' '<<g->linMom[0]<<' '<<g->linMom[1]<<' '<<g->angMom<<' '<<g->eKE<<' '<<g->eBForce<<' '<<g->eElast<<' '<<g->eElastExact<<' '<<g->eViscF<<' '<<g->eViscB<<' '<<g->eFF<<' '<<g->eTotal<<' '<<g->wtTotalOutStep.tv_sec+g->wtTotalOutStep.tv_usec/1000000.0<<' '<<g->wtTotalForceCalc.tv_sec+g->wtTotalForceCalc.tv_usec/1000000.0<<' '<<g->wtTotalDddtCalc.tv_sec+g->wtTotalDddtCalc.tv_usec/1000000.0<<' '<<g->wtTotalReset.tv_sec+g->wtTotalReset.tv_usec/1000000.0<<' '<<g->wtTotalUpdateDomain.tv_sec+g->wtTotalUpdateDomain.tv_usec/1000000.0<<' '<<g->wtTotalMPI.tv_sec+g->wtTotalMPI.tv_usec/1000000.0<<' '<<g->maxV<<' '<<g->aveDens<<' '<<g->varDens<<' '<<g->aveDensFromMass<<' '<<g->n<<' '<<g->nSph<<' '<<g->maxF<<' '<<g->maxFF<<' '<<g->rmsFF<<' '<<g->edElastdt<<' '<<g->edViscFdt<<' '<<g->edFFdt;
   for (int i=0;i<GLOBAL_CUSTOM_BUFFER_SIZE;i++) {
      fio <<' '<<g->custom[i];
   }
   fio << '\n';
}

void CioGlobals::openAndSeekToStep(int outStep,CglobalVars *g) {
   fio.open(filename.c_str(),ios::in|ios::out);

   if (!fio.is_open() ) {
      cout <<"Error: globals file will not open!!!!"<<endl;
      exit(-1);
   }

   string dummy;
   getline(fio,dummy);
   fio.seekp(fio.tellg());
   fio.clear();

   if (outStep > 0) {
      double currStep = myGetLine(g);
      while ((currStep < outStep)&&!fio.eof()) {
         currStep = myGetLine(g);
      }
      if (currStep != outStep) {
         fio.seekp(0,ios::end);
         fio.clear();
      } else {
         fio.seekp(fio.tellg());
      }
   }
}


int CioGlobals::findClosestStepToTime(double time) {
   fstream f;
   f.open(filename.c_str(),ios::in);
   if (!f.is_open() ) {
      cout <<"Error: globals file will not open!!!!"<<endl;
      exit(-1);
   }
   string line;
   getline(f,line);
   getline(f,line);
   istringstream ss(line);
   int currStep,sphStep;
   double thisTime,lastTime;

   ss >> currStep >> sphStep >> thisTime;
   lastTime = thisTime;
   while (thisTime<time) {
      lastTime = thisTime;
      getline(f,line);
      ss.str(line);
      ss >> currStep >> sphStep >> thisTime;
   }
   if (abs(lastTime-time)<abs(thisTime-time)) {
      return currStep-1;
   } else {
      return currStep;
   }
}


void CioGlobals::readGlobals(int outStep,CglobalVars *g) {
   openAndSeekToStep(outStep,g);
   fio.close();
}       

void CioGlobals::writeGlobals(int outStep,CglobalVars *g) {
   CglobalVars sendBuf = *g;
   vector<CglobalVars> recvBuf;
   int gSize = sizeof(CglobalVars);
   if (g->mpiRank==0) {
      recvBuf.resize(g->mpiSize);
   }

   if (g->mpiSize > 1) {
      MPI_Gather(&sendBuf,gSize,MPI_BYTE,&(recvBuf[0]),gSize,MPI_BYTE,0,MPI_COMM_WORLD);
   } else if (g->mpiRank==0) {
      recvBuf[0] = *g;
   }

   CglobalVars tmpg;
   openAndSeekToStep(outStep-1,&tmpg);
   if (g->mpiRank==0) {
      CglobalVars sumG;
      sumG.init(recvBuf[0]);
      for (int i=0;i<g->mpiSize;i++) {
         sumG += recvBuf[i];
      }
      cout << "\tWriting Global Vars with outStep = "<<outStep <<endl;
      myWriteLine(outStep,&sumG);
      recvBuf.clear();
   } else {
      myWriteLine(outStep,g);
   }
   fio.close();
}

