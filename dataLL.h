
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


class CpInfoLink;

/*
 * CpInfo contains information on a particle (Cparticle) that is used by
 * CdataLL. This information is: 
 *            - key        used to spatially sort the particles
 *            - neighbrs   list of neighbours     
 *            - ppInfoLink link to entry for particle in pInfoLink.
 */
class CpInfo {
public:
      Cparticle *p;
      CpInfoLink *ppInfoLink;
      unsigned long key;
      vector<Cparticle *> neighbrs;


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

/*
 * contains back-link to particle information. Used to implement vector of
 * pointers pInfoLinks that exists alongside ps, the vector of particles. pInfoLinks points to particle
 * information in pInfos, a list that is spatially sorted. This way the CdataLL
 * class can loop through a sorted list of particle (pInfos) or an unsorted
 * list (ps or pInfoLinks). 
 */
class CpInfoLink {
public:
      CpInfo *ppInfo;
};

/*
* helper classes for the functions
*     CdataLL::getParticleFromTag(int tag);
*     CdataLL::getPinfoFromTag(int tag);
*/
class CtagPP {
public:
     int tag;
     Cparticle *p;
     CtagPP() {
        tag = -1;
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


/*
 * the main particle data class. This handles:
 *    -storing the particle data in ps
 *    -applying functions to particles (traverse())
 *    -applying functions to particles and their neighbours (neighbours() and neighboursGroups())
 *    -providing a sorted (pInfos) and unsorted (pInfoLinks and ps) list of
 *        particles and particle information (neighbours etC). 
 *        Sorting done spatially to improve cache performance
 *    -mpi communication of particle and ghost particles to neighbouring cpus
 *    -get particle by tag
 *    -misc others....
 */
class CdataLL {
public:

      /*
       * constructor takes a vector of particles (in_ps), global variables
       * (in_globals) and a flag that determines whether spatial ordering is to be performed
       */
      CdataLL(particleContainer &in_ps, CglobalVars &in_globals,bool opt);
      
      /*
       * reset - Whenever the particle positions are changed the neighbour
       * lists and the ghost particles from neighbouring cpus are no longer valid.
       * This function resets these, and all other data contained in dataLL
       * that depends on the particle positions.
       * The user would normally call this after the particles have been moved.
       */
      void reset();

      /*
       * traverse - This applies a function (Thefunct) to all the particles.
       *            [optional] -can supply a function (Iffunct) that only
       *                        applies Thefunct if Iffunct returns true.
       *                       -can supply input constants or state variables.
       *                        see below for details....
       *
       * big TODO: all these different options have evolved over the life of
       *           this code. They can, and should, be combined into a smaller
       *           number (optimally just one) of functions. 
       */
      template <void Thefunct(Cparticle &,CglobalVars &)>
      void traverse();
      template <void Thefunct(Cparticle &,CglobalVars &,double)>
      void traverse(double);
      template <class T, void Thefunct(Cparticle &,CglobalVars &,T &)>
      void traverse(T &);
      template <void Thefunct(Cparticle &,CglobalVars &), bool Iffunct(Cparticle &)>
      void traverse();
      template <void Thefunct(Cparticle &,CglobalVars &,double), bool Iffunct(Cparticle &)>
      void traverse(double);
      template <class T, void Thefunct(Cparticle &,CglobalVars &,T &), bool Iffunct(Cparticle &)>
      void traverse(T &);
      template <void Thefunct(Cparticle &,CglobalVars &,double,double)>
      void traverse(double,double);
      template <void Thefunct(Cparticle &,CglobalVars &,double,double), bool Iffunct(Cparticle &)>
      void traverse(double,double);
      template <void Thefunct(Cparticle &,CglobalVars &,double,double,double), bool Iffunct(Cparticle &)>
      void traverse(double,double,double);


      /*
       * neighbours - This applies a function (Thefunct) to all the particle
       *              pairs that are neighbours with each other
       *              [optional] -can supply a function (Iffunct) that only
       *                        applies Thefunct if Iffunct returns true.
       */
      template <void Thefunct(Cparticle &,Cparticle &,CglobalVars &)>
      void neighbours();
      template <void Thefunct(Cparticle &,Cparticle &,CglobalVars &), bool Iffunct(Cparticle &)>
      void neighbours();

      /*
       * neighboursGroup - This applies a function (Thefunct) to all the neighbours
       *              of each particle. The neighbours are passed to the
       *              function using a vector of particle pointers.
       *              [optional] -can supply a function (Iffunct) that only
       *                          applies Thefunct if Iffunct returns true.
       *                         -can supply a state variable (with an
       *                          arbitrary class)
       */
      template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &)>
      void neighboursGroup();
      template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &), bool Iffunct(Cparticle &)>
      void neighboursGroup();
      template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &), bool Iffunct(Cparticle &)>
      void neighboursGroupAtRadius(const double radius);
#ifdef LIQ_DEM
      template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &), bool Iffunct(Cparticle &)>
      void neighboursGroupCoupling();
#endif
      template <class T, void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &,T &), bool Iffunct(Cparticle &)>
      void neighboursGroup(T &);

      /*
       * neighboursUsing - This applies a function (Thefunct) to all the particles
       *              that are neighbours with each particle in _ps.
       *              [optional] -Thefunct can either take a particle pair or a
       *                          particle and all its neighbours (as a vector of particle pointers).
       */
      template <void Thefunct(Cparticle &,Cparticle &,CglobalVars &)>
      void neighboursUsing(vector<Cparticle> &_ps);
      template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &)>
      void neighboursUsing(vector<Cparticle> &_ps);
      
      /*
       * syncParticlesBetweenProcs - syncs data between neighbouring cpus. The
       * data to be synced is read/writen using the templated functions (to be
       * applied to each particles).
       */
      template <class T, T readFunct(Cparticle &), void writeFunct(Cparticle &,T)>
      void syncParticlesBetweenProcs();
      template <class T, T readFunct(Cparticle &), void writeFunct(Cparticle &,T)>
      void reverseSyncParticlesBetweenProcs();

      /*
       * functionOverGrid - apply theFunct over a set of grid points and the
       * particles that neighbour these grid points. theFunct
       * takes as input a particle pair, the first particle is a grid point, 
       * the second is a neighbour of that grid point. The final grid points
       * are output in outPs
       */
      template <void theFunct(Cparticle &,Cparticle &,CglobalVars &)>
      void functionOverGrid(vector<Cparticle> &outPs,vectInt &gridDimN);

      /*
       * find neighbours of particle _p. neighbours returned in _neighbrs
       */
      void calcNeighbours(vector<Cparticle *> &_neighbrs,Cparticle &_p);
      void calcNeighboursAtRadius(vector<Cparticle *> &_neighbrs,Cparticle &_p,const double radius,const double radius2);
#ifdef LIQ_DEM
      void calcNeighboursCoupling(vector<Cparticle *> &_neighbrs,Cparticle &_p);
#endif


      /*
       * set a timestep (sent to all cpus, the minimum timestep given by all
       * cpus is used)
       */
      void setGlobalTimestep(double dt);

      /*
       * print info on dataLL class
       */
      void printMe();

      /*
       * update the spatial ordering of the particles. This improves cache
       * performance. Only use once particles are sufficiently disordered.
       * (eg. every output timestep)
       */
      void calcTraverseOrder();

      int numParticles() { return ps.size(); };
      particleContainer *getParticles() { return &ps; }

      /*
       * useful if you have just made massive changes to the particle array 
       * in general or particle positions in general. Reconstructs the spatial
       * ordering from scratch (warning: very slow!)
       */
      void updateParticles();

      /*
       * insertion/deletion functions
       */
      void insertNewParticle(Cparticle &newP);
      void markForDeletion(Cparticle &p);
      void deleteParticles();

      /*
       * get particle or particle info from a (unique) tag
       */
      Cparticle *getParticleFromTag(int tag);
      CpInfo *getPinfoFromTag(int tag);

      /*
       * sum a value over all cpus
       */
      void sumOverProcs(vect *data,int num);
      void sumOverProcs(double *data,int num);

      /*
       * global variables. TODO: this should be private
       */
      CglobalVars &globals;

private:

      /*
       * functions to allocate mpi buffers
       */
      int calcBufferSize(vectInt coords);
      int calcPBufferSize(vectInt coords);
      void allocateBuffers();

      /*
       * delete particle from ps, pInfos and pInfoLink
       */
      void deleteParticle(vector<CpInfoLink>::iterator &ppInfoLink,CpInfo *ppInfo);

      /*
       * insert a particle into the bucket array
       */
      void insert(Cparticle &);

      void calcDomainLimits();

      /*
       * add _toAdd to list of neighbours of _p, only if _toAdd is a neighbour
       * of _p
       */
      void addIfNeighbr(Cparticle &_p,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd);
      void addIfInRadius(Cparticle &_p,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd,const double radius2);
      void addAll(Cparticle &_p,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd);

      /*
       * maps from position to bucket cell index
       */
      inline vectInt getCellI(vect r) {
         return static_cast<vectInt>((r-cell_min)/(KERNAL_RADIUS*sph_search_radius)+3);
      }

#ifdef LIQ_DEM
      inline vectInt dem_getCellI(vect r) {
         //return static_cast<vectInt>((r-rmin)/(2*DEM_SEARCH_RADIUS)+3*H/DEM_SEARCH_RADIUS-1);
         return static_cast<vectInt>((r-cell_min)/(2*DEM_SEARCH_RADIUS)+2*LIQ_DEM_COUPLING_RADIUS/DEM_SEARCH_RADIUS-1);
      }
#endif

      /*
       * functions to sort particles spatially
       *    calcKey. calculates a key for each particle using morton ordering
       *    initPinfos - clears current pInfos (vector of particle information,
       *    like neighbours, keys etc. and recreates it
       *    checkPinfos - checks that pInfos is consistent with pInfoLinks
       *    sortPInfos - sorts pInfos according to key
       *    initTraverseOrder - runs previous three functions.
       */
      unsigned long calcKey(vect r);
      void initPinfos();
      void checkPinfos();
      void sortPInfos();
      void initTraverseOrder(bool opt);
 
      /*
       * functions that implement mpi communication
       */
      void updateDomain();
      void sendRecvGhosts(const vectInt coords,const int num,MPI_Request *request,int *ghostSizesSend,int *ghostSizesRecv);
      void sendRecvParticles(const vectInt coords,const int num,MPI_Request *request,int *particleSizesSend,int *particleSizesRecv);
      void sendRecvDataSync(const vectInt coords,const int num,MPI_Request *request,int sendSize, void *sendBuffer, int recvSize, void *recvBuffer);
      
      /*
       * vectors that hold pointers to particles or pInfos. These pointers are
       * sorted by tag, so that the particle or pInfo with a given tag can easily be found
       */
      vector<CtagPP> tagSortedParticlePointers;
      vector<CtagPI> tagSortedPinfoPointers;


      /*
       * domain extents and max smoothing length
       */
      vect rmax;
      vect rmin;
      double hmax;
      double sph_search_radius;

      Array<vector<Cparticle *>,NDIM> cells; //bucket list
      vect cell_min,cell_max;
      vector<vector<Cparticle *>*> dirty_cells;
#ifdef LIQ_DEM
      Array<vector<Cparticle *>,NDIM> dem_cells; //bucket list
#endif
      int n;                                 //number of particles
      particleContainer &ps;                 //holds the set of particles
      Array<Cparticle **,NDIM> ghostedParticles; //points to the particles ghosted to neighbouring cpus
      vector<Cparticle> procGhostParticles;  //ghost particles received from neighbouring cpus
      Array<double,NDIM> procGhostMax2H;     //max 2h from neighbouring cpus


      /*
       * buffers and buffer sizes for mpi communication
       */
      Array<CghostData*,NDIM> ghostBuffersSend;
      Array<CghostData*,NDIM> ghostBuffersRecv;
      Array<Cparticle*,NDIM> particleBuffersSend;
      Array<Cparticle*,NDIM> particleBuffersRecv;
      Array<void*,NDIM> syncBuffersSend;
      Array<void*,NDIM> syncBuffersRecv;

      Array<int,NDIM> pBufferSizes;
      Array<int,NDIM> gBufferSizes;
      Array<int,NDIM> sendSizesGhosts;
      Array<int,NDIM> recvSizesGhosts;

      /*
       * pInfos holds particle information that is used by this class
       * (neighbours, sorting key etc.). 
       * pInfoLinks is a vector where the i'th entry contains a pointer to the
       * pInfo for particle i in the main particle vector (ps)
       */
      vector<CpInfoLink> pInfoLinks;
      vector<CpInfo> pInfos;
};

#endif
