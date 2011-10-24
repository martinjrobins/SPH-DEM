
#ifndef DATALL_H
#define DATALL_H

#include <vector>
#include <list>
#include <blitz/array.h>
#include <mpi.h>
#include <sys/time.h>
#include "customConstants.h"
#include "particle.h"
#include "sphConstants.h"
#include "globalVars.h"




class CdataLL {
public:
      int calcBufferSize(vectInt coords);
      CdataLL(vector<particleT> &in_ps, CglobalVars &in_globals,bool opt): ps(in_ps),globals(in_globals) {
         hmax = H;
         gBufferSize = MAX_NUM_PARTICLES_PER_CPU;
         pBufferSize = MAX_NUM_PARTICLES_PER_CPU;
         //gBufferSize = int(ps.size()*4/(pow(3.0,NDIM)-2));
         //pBufferSize = int(ps.size()*4/(10*pow(3.0,NDIM)));
         cout << "gBufferSize = "<<gBufferSize<<" pBufferSize = "<<pBufferSize<<endl;
         ghostBuffersSend.resize(3);

         for (Array<CghostData*,NDIM>::iterator i=ghostBuffersSend.begin();i!=ghostBuffersSend.end();i++) {
            vectInt coords = i.position()-1;
            if (globals.procNeighbrs(i.position()) >= 0) {
               *i = new CghostData[calcBufferSize(coords)];
            }
         }
         particleBuffersSend.resize(3);
         for (Array<Cparticle*,NDIM>::iterator i=particleBuffersSend.begin();i!=particleBuffersSend.end();i++) {
            vectInt coords = i.position()-1;
            if (globals.procNeighbrs(i.position()) >= 0) {
               *i = new Cparticle[calcBufferSize(coords)/4];
            }
         }
         ghostBuffersRecv.resize(3);
         for (Array<CghostData*,NDIM>::iterator i=ghostBuffersRecv.begin();i!=ghostBuffersRecv.end();i++) {
            vectInt coords = i.position()-1;
            if (globals.procNeighbrs(i.position()) >= 0) {
               *i = new CghostData[calcBufferSize(coords)];
            }
         }
         particleBuffersRecv.resize(3);
         for (Array<Cparticle*,NDIM>::iterator i=particleBuffersRecv.begin();i!=particleBuffersRecv.end();i++) {
            vectInt coords = i.position()-1;
            if (globals.procNeighbrs(i.position()) >= 0) {
               *i = new Cparticle[calcBufferSize(coords)/4];
            }
         }
         ghostedParticles.resize(3);
         for (Array<Cparticle**,NDIM>::iterator i=ghostedParticles.begin();i!=ghostedParticles.end();i++) {
            vectInt coords = i.position()-1;
            if (globals.procNeighbrs(i.position()) >= 0) {
               *i = new Cparticle*[calcBufferSize(coords)];
            }
         }
         syncBuffersSend.resize(3);
         for (Array<void*,NDIM>::iterator i=syncBuffersSend.begin();i!=syncBuffersSend.end();i++) {
            vectInt coords = i.position()-1;
            if (globals.procNeighbrs(i.position()) >= 0) {
               *i = new unsigned char[sizeof(CghostData)*calcBufferSize(coords)];
            }
         }
         syncBuffersRecv.resize(3);
         for (Array<void*,NDIM>::iterator i=syncBuffersRecv.begin();i!=syncBuffersRecv.end();i++) {
            vectInt coords = i.position()-1;
            if (globals.procNeighbrs(i.position()) >= 0) {
               *i = new unsigned char[sizeof(CghostData)*calcBufferSize(coords)];
            }
         }

         sendSizesGhosts.resize(3);
         recvSizesGhosts.resize(3);
         procGhostMax2H.resize(3);
         procGhostMax2H = 1.1*KERNAL_RADIUS*H;

         cout << "Processor "<<globals.mpiRank<<" of "<<globals.mpiSize<<" has domain ";
         for (int i=0;i<NDIM*2;i++) {
            cout<<globals.procDomain[i]<<' ';
         }

#ifdef GRID_SMOOTHING
         vectInt smSize;
         for (int i=0;i<NDIM;i++) {
            double size = (RMAX[i]-RMIN[i])/(SMOOTH_GRID_SIZE);
            smSize[i] = int(floor((RMAX[i]-RMIN[i])/(SMOOTH_GRID_SIZE) + 0.5));
            if (smSize[i]/size != 1) {
               cerr << "Error: SMOOTH_GRID_SIZE doesn't give an integer number of grid cells in domain!" << endl;
               exit(-1);
            }
         }  
         vectInt minI,maxI;
         for (int i=0;i<NDIM;i++) {
            double test = (globals.procDomain[i*2]-RMIN[i])/SMOOTH_GRID_SIZE;
            minI[i] = int(floor(test));
            if (test==double(minI[i])) {
               minI[i]--;
            }
            minI[i]--;
            test = (globals.procDomain[i*2+1]-RMIN[i])/SMOOTH_GRID_SIZE;
            maxI[i] = int(ceil(test));
            maxI[i]++;
            smSize[i] = maxI[i]-minI[i]+1;
            smGridMin[i] = minI[i]*SMOOTH_GRID_SIZE+RMIN[i];
         }
         cout <<"Processor "<<globals.mpiRank<<": smSize = "<<smSize<<"smGridMin = "<<smGridMin<<" minI = "<<minI<<" maxI = "<<maxI<<endl;

         smGrid.resize(smSize);
         for (Array<CsmVertex,NDIM>::iterator i=smGrid.begin();i!=smGrid.end();i++) {
            vectInt coords = i.position();
            for (int j=0;j<NDIM;j++) {
               i->pos[j] = coords[j]*SMOOTH_GRID_SIZE + smGridMin[j];
            }
         }

         smGridBuffersRecv.resize(3);
         smGridBuffersSend.resize(3);
         smGridViews.resize(3);
         for (Array<CsmVertex*,NDIM>::iterator i = smGridBuffersSend.begin();i != smGridBuffersSend.end();i++) {
            vectInt coords = i.position();
            int neighbr = globals.procNeighbrs(coords);
               if (neighbr >= 0) {
               vectInt upperBounds,lowerBounds;
               for (int j=0;j<NDIM;j++) {
                  int extent = smGrid.extent(j);
                  if (coords[j]==1) {
                     lowerBounds[j] = 0;
                     upperBounds[j] = extent-1;
                  } else if (coords[j]==0) {
                     lowerBounds[j] = 0;
                     upperBounds[j] = 3;
                  } else if (coords[j]==2) {
                     lowerBounds[j] = extent-4;
                     upperBounds[j] = extent-1;
                  }
               
               }
               RectDomain<NDIM> subdomain(lowerBounds, upperBounds);
               smGridViews(coords).reference(smGrid(subdomain));
               *i = new CsmVertex[smGridViews(coords).size()];
               smGridBuffersRecv(coords) = new CsmVertex[smGridViews(coords).size()];
            }
         }
#endif
         if (opt) {
            initTraverseOrder();
         }
         reset();
      }
      void reset();
      void reset(const double scale);

      void sumOverProcs(vect *data,int num);
      void sumOverProcs(double *data,int num);
      template <class T, T readFunct(particleT &), void writeFunct(particleT &,T)>
      void syncParticlesBetweenProcs();
      template <void Thefunct(particleT &,CglobalVars &)>
      void traverse();
      template <class T, void Thefunct(particleT &,CglobalVars &,T &)>
      void traverse(T &);
      template <void Thefunct(particleT &,CglobalVars &), bool Iffunct(particleT &)>
      void traverse();
      template <class T, void Thefunct(particleT &,CglobalVars &,T &), bool Iffunct(particleT &)>
      void traverse(T &);
      template <void Thefunct(particleT &,particleT &,CglobalVars &)>
      void neighboursUsing(vector<particleT> &_ps);
      template <void Thefunct(particleT &,vector<particleT *> &,CglobalVars &)>
      void neighboursUsing(vector<particleT> &_ps);
      template <void Thefunct(particleT &,particleT &,CglobalVars &)>
      void neighbours();
      template <void Thefunct(particleT &,particleT &,CglobalVars &), bool Iffunct(particleT &)>
      void neighbours();
      template <void Thefunct(particleT &,vector<particleT *> &,CglobalVars &)>
      void neighboursGroup();
      template <void Thefunct(particleT &,vector<particleT *> &,CglobalVars &), bool Iffunct(particleT &)>
      void neighboursGroup();
      template <class T, void Thefunct(particleT &,vector<particleT *> &,CglobalVars &,T &), bool Iffunct(particleT &)>
      void neighboursGroup(T &);
      template <void theFunct(particleT &,particleT &,CglobalVars &)>
      void functionOverGrid(vector<particleT> &outPs,vectInt &gridDimN);
      void calcNeighbours(vector<particleT *> &_neighbrs,particleT &_p);

      void setGlobalTimestep(double dt);
      void printMe();
      void calcTraverseOrder();
      int numParticles() { return ps.size(); };
      vector<particleT> *getParticles() { return &ps; }
      void updateParticles();




      void insertNewParticle(particleT &newP);

      void markForDeletion(particleT &p);
      void deleteParticles();

      particleT *getParticleFromTag(int tag);
      CpInfo *getPinfoFromTag(int tag);

      template <void Thefunct(particleT &,vector<particleT *> &,CglobalVars &), bool Iffunct(particleT &)>
      void neighboursGroupScaled();


      template <bool Iffunct(particleT &)>
      void findBothNeighbours();

      template <bool Iffunct(particleT &)>
      void findNeighbours();

      template<void Thefunct(particleT &,Array<CsmVertex*,NDIM> &,CglobalVars &), bool Iffunct(particleT &)>
      //template<void Thefunct(particleT &,vector<CsmVertex *> &,CglobalVars &), bool Iffunct(particleT &)>
      void useSmGrid();
      template<void Thefunct(CsmVertex &,CglobalVars &)>
      void traverseSmGrid();

      void syncSmGrid();

      CglobalVars &globals;
private:

      class CsmVertex {
      public:
         vect v;
         vect pos;
         int num;
      };

      class CtagPP {
      public:
         int tag;
         particleT *p;
         CtagPP() {
            tag = -1;
         }
      };

      class CpInfoLink;

      class CpInfo {
      public:
         particleT *p;
         CpInfoLink *ppInfoLink;
         unsigned long key;
         vector<particleT *> neighbrs;
         vector<particleT *> scaledNeighbrs;

         CpInfo() {}

         CpInfo(const CpInfo &pInfo) {
            p = pInfo.p;
            ppInfoLink = pInfo.ppInfoLink;
            key = pInfo.key;
         }

         CpInfo& operator=(CpInfo &pInfo) {
            p = pInfo.p;
            ppInfoLink = pInfo.ppInfoLink;
            key = pInfo.key;
            return *this;
         }
      };

      class CtagPI {
      public:
            int tag;
            CpInfo *pInfo;
            CtagPI() {
               tag = -1;
            }
      };

      class CpInfoLink {
      public:
            CpInfo *ppInfo;
      };


      void insertionSort(std::list<CpInfo> &v);

      void sendRecvSmGridData(vectInt coords,int num,MPI_Request *requestSendRecv,CsmVertex *sendBuf,CsmVertex *recvBuf,int size);
      void getVerticesFromPosition(vect pos,vectInt &lowerB);
      void calcDomainLimits();
      void insert(particleT &,const double);
      void addIfNeighbr(particleT &_p,vector<particleT *> &_neighbrs,vector<particleT *> &_toAdd);
      void deleteParticle(vector<CpInfoLink>::iterator &ppInfoLink,CpInfo *ppInfo);

      void addIfBothNeighbr(particleT &_p,vector<particleT *> &_scaledNeighbrs,vector<particleT *> &_neighbrs,vector<Cparticle *> &_toAdd);
      void calcBothNeighbours(vector<particleT *> &_scaledNeighbrs,vector<particleT *> &_neighbrs,particleT &_p);
      void addIfScaledNeighbr(particleT &_p,vector<particleT *> &_neighbrs,vector<particleT *> &_toAdd);
      void calcScaledNeighbours(vector<particleT *> &_neighbrs,particleT &_p);

      inline vectInt getCellI(vect r) {
         return static_cast<vectInt>((r-rmin)/(KERNAL_RADIUS*scale*hmax)+2);
      }

      void insertGhosts();
      void removeGhosts();
      unsigned long calcKey(vect r);
      void initPinfos();
      void checkPinfos();
 
      void initTraverseOrder();
      void updateDomain();
      void sendRecvData(const vectInt coords,const int num,MPI_Request *request,int *ghostSizesSend,int *ghostSizesRecv,int *particleSizesSend,int *particleSizesRecv);
      void sendRecvGhosts(const vectInt coords,const int num,MPI_Request *request,int *ghostSizesSend,int *ghostSizesRecv);
      void sendRecvParticles(const vectInt coords,const int num,MPI_Request *request,int *particleSizesSend,int *particleSizesRecv);
      void sendRecvDataSync(const vectInt coords,const int num,MPI_Request *request,int sendSize, void *sendBuffer, int recvSize, void *recvBuffer);
      void sortPInfos();
      
      vector<CtagPP<particleT> > tagSortedParticlePointers;
      vector<CtagPI<particleT> > tagSortedPinfoPointers;
      double scale;
      vector<particleT> &ps;
      vect rmax;
      vect rmin;
      double hmax;
      //vector<vector<particleT *> > neighbrs;
      Array<vector<particleT *>,NDIM> cells;
      int n;
      vector<CpInfo> pInfos;
      //vector<particleT *> traverseOrder;
      //vector<unsigned long> traverseKeys;
      Array<particleT **,NDIM> ghostedParticles;
      //vector<particleT>::iterator startOfRecvGhostParticles;

      Array<void*,NDIM> syncBuffersSend;
      Array<void*,NDIM> syncBuffersRecv;
      Array<int,NDIM> sendSizesGhosts;
      Array<int,NDIM> recvSizesGhosts;

      Array<particleT::ghost*,NDIM> ghostBuffersSend;
      Array<particleT*,NDIM> particleBuffersSend;
      Array<CsmVertex*,NDIM> smGridBuffersSend;
      Array<particleT::ghost*,NDIM> ghostBuffersRecv;
      Array<particleT*,NDIM> particleBuffersRecv;
      Array<CsmVertex*,NDIM> smGridBuffersRecv;
      Array<Array<CsmVertex,NDIM>,NDIM> smGridViews;
      int pBufferSize;
      int gBufferSize;

      vector<particleT> procGhostParticles;
      Array<double,NDIM> procGhostMax2H;
      Array<CsmVertex,NDIM> smGrid;
      vect smGridMin;

      vector<CpInfoLink> pInfoLinks;
};

#include "dataLL.impl"

#endif
