#include "ioDomain.h"

void CioDomain::readDomain(int outStep,CglobalVars *g) {
   //open file for reading
   fio.open(filename.c_str(),ios::in);
   char buffer[200];
   fio >> g->mpiRank;
   for (int i=0;i<NDIM*2;i++) {
      fio >> g->procDomain[i];
   }
   for (Array<int,NDIM>::iterator p=g->procNeighbrs.begin();p!=g->procNeighbrs.end();p++) {
      fio >> *p;
   }
   vectInt tmp = 1;
   g->procNeighbrs(tmp) = -1;
   
   cout << "Processor "<<g->mpiRank<<" of "<<g->mpiSize<<" has domain ";
   for (int i=0;i<NDIM*2;i++) {
      cout<<g->procDomain[i]<<' ';
   }
   cout << " and neighbrs "<<g->procNeighbrs<<endl;
   fio.close();
}


void CioDomain::writeDomain(int outStep,CglobalVars *g) {
   //open file (overwrite existing contents)
   fio.open(filename.c_str(),ios::out|ios::trunc);

   //write domain
   fio << g->mpiRank << ' ';
   for (int j=0;j<NDIM*2;j++) {
      fio << g->procDomain[j] << ' ';
   }
   for (Array<int,NDIM>::iterator p=g->procNeighbrs.begin();p!=g->procNeighbrs.end();p++) {
      fio << *p << ' ';
   }
   fio << endl;
   fio.close();
}


void CioDomain::readDomainOld(int outStep,CglobalVars *g) {
   //open file for reading
   fio.open(filename.c_str(),ios::in);
   char buffer[200];
   for (int i=0;i<g->mpiRank;i++) {
      fio.getline(buffer,200); 
   }
   int ckRank;
   fio >> ckRank;

   for (int i=0;i<NDIM*2;i++) {
      fio >> g->procDomain[i];
   }
   for (Array<int,NDIM>::iterator p=g->procNeighbrs.begin();p!=g->procNeighbrs.end();p++) {
      fio >> *p;
   }
   vectInt tmp = 1;
   g->procNeighbrs(tmp) = -1;
   
   cout << "Processor "<<g->mpiRank<<" of "<<g->mpiSize<<" has domain ";
   for (int i=0;i<NDIM*2;i++) {
      cout<<g->procDomain[i]<<' ';
   }
   cout << " and neighbrs "<<g->procNeighbrs<<endl;
    

   fio.close();
}

void CioDomain::writeDomainOld(int outStep,vector<vector<double> > vprocDomain, vector<Array<int,NDIM> > vprocNeighbrs) {
   //open file (overwrite existing contents)
   fio.open(filename.c_str(),ios::out|ios::trunc);

   //write domain
   int n = vprocDomain.size();
   for (int i=0;i<n;i++) {
      int m = vprocDomain[i].size();
      fio << i << ' ';
      for (int j=0;j<m;j++) {
         fio << vprocDomain[i][j] << ' ';
      }

      for (Array<int,NDIM>::iterator p=vprocNeighbrs[i].begin();p!=vprocNeighbrs[i].end();p++) {
         fio << *p << ' ';
      }
      
      fio << endl;
   }
   fio.close();
}


void CioDomain::setFilename(string _filename,CglobalVars *g) {
   char str[20];
   sprintf(str,"Domain%d.dat",g->mpiRank);
   filename = _filename+str;
}
