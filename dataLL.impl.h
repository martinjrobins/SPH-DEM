#ifndef DATALL_IMPL_H
#define DATALL_IMPL_H

#include "dataLL.h"
#include "misc.h"


/*
 * basic traverse function. Use this to apply a per-particle function to the
 * set of particles (eg. a drift or kick operation)
 * loops through the ps vector (vector of Cparticle) and runs T_thefunct on
 * each particle
 * also loops through all the ghost particles (procGhostParticles) and does
 * the same.
 */
template <void T_thefunct(Cparticle &,CglobalVars &)>
void CdataLL::traverse() {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      T_thefunct(*pp,globals);
   }
   for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
      T_thefunct(*pp,globals);
   }
}

/*
 * same as above but adds an input constant to T_thefunct
 */
template <void T_thefunct(Cparticle &,CglobalVars &, double)>
void CdataLL::traverse(double param) {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      T_thefunct(*pp,globals,param);
   }
   for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
      T_thefunct(*pp,globals,param);
   }
}

/*
 * same as above but only applies the function if the boolean Iffunct is true
 */
template <void T_thefunct(Cparticle &,CglobalVars &, double),bool Iffunct(Cparticle &)>
void CdataLL::traverse(double param) {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      if (Iffunct(*pp)) {
         T_thefunct(*pp,globals,param);
      }
   }
   //TODO: ghost or not should be indicated by another flag, but this would be a significant change
   if (!procGhostParticles.empty()) {
      enum iamTypes store = procGhostParticles.begin()->iam;
      procGhostParticles.begin()->iam = ghost;
      if (Iffunct(*(procGhostParticles.begin()))) {
         procGhostParticles.begin()->iam = store;
         for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
            if (Iffunct(*pp)) {
               T_thefunct(*pp,globals,param);
            }
         }
      }
      procGhostParticles.begin()->iam = store;
   }
}

/*
 * and other variants. Look at header file for info...
 */
template <class T, void Thefunct(Cparticle &,CglobalVars &,T &)>
void CdataLL::traverse(T &param) {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      Thefunct(*pp,globals,param);
   }
   for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
      Thefunct(*pp,globals,param);
   }
}

template <class T, void Thefunct(Cparticle &,CglobalVars &,T &), bool Iffunct(Cparticle &)>
void CdataLL::traverse(T &param) {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      if (Iffunct(*pp)) {
         Thefunct(*pp,globals,param);
      }
   }
   if (!procGhostParticles.empty()) {
      enum iamTypes store = procGhostParticles.begin()->iam;
      procGhostParticles.begin()->iam = ghost;
      if (Iffunct(*(procGhostParticles.begin()))) {
         procGhostParticles.begin()->iam = store;
         for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
            if (Iffunct(*pp)) {
               Thefunct(*pp,globals,param);
            }
         }
      }
      procGhostParticles.begin()->iam = store;
   }
}

template <void T_thefunct(Cparticle &,CglobalVars &, double, double)>
void CdataLL::traverse(double param1,double param2) {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      T_thefunct(*pp,globals,param1,param2);
   }
   for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
      T_thefunct(*pp,globals,param1,param2);
   }
}

template <void T_thefunct(Cparticle &,CglobalVars &, double, double),bool Iffunct(Cparticle &)>
void CdataLL::traverse(double param1, double param2) {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      if (Iffunct(*pp)) {
         T_thefunct(*pp,globals,param1,param2);
      }
   }
   if (!procGhostParticles.empty()) {
      enum iamTypes store = procGhostParticles.begin()->iam;
      procGhostParticles.begin()->iam = ghost;
      if (Iffunct(*(procGhostParticles.begin()))) {
         procGhostParticles.begin()->iam = store;
         for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
            if (Iffunct(*pp)) {
               T_thefunct(*pp,globals,param1,param2);
            }
         }
      }
      procGhostParticles.begin()->iam = store;
   }
}

template <void T_thefunct(Cparticle &,CglobalVars &, double, double, double),bool Iffunct(Cparticle &)>
void CdataLL::traverse(double param1, double param2, double param3) {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      if (Iffunct(*pp)) {
         T_thefunct(*pp,globals,param1,param2,param3);
      }
   }
   if (!procGhostParticles.empty()) {
      enum iamTypes store = procGhostParticles.begin()->iam;
      procGhostParticles.begin()->iam = ghost;
      if (Iffunct(*(procGhostParticles.begin()))) {
         procGhostParticles.begin()->iam = store;
         for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
            if (Iffunct(*pp)) {
               T_thefunct(*pp,globals,param1,param2,param3);
            }
         }
      }
      procGhostParticles.begin()->iam = store;
   }
}

template <void T_thefunct(Cparticle &,CglobalVars &),bool Iffunct(Cparticle &)>
void CdataLL::traverse() {
   for (particleContainer::iterator pp = ps.begin();pp!=ps.end();pp++) {
      if (Iffunct(*pp)) {
         T_thefunct(*pp,globals);
      }
   }
   if (!procGhostParticles.empty()) {
      enum iamTypes store = procGhostParticles.begin()->iam;
      procGhostParticles.begin()->iam = ghost;
      if (Iffunct(*(procGhostParticles.begin()))) {
         procGhostParticles.begin()->iam = store;
         for (vector<Cparticle>::iterator pp = procGhostParticles.begin();pp!=procGhostParticles.end();pp++) {
            if (Iffunct(*pp)) {
               T_thefunct(*pp,globals);
            }
         }
      }
      procGhostParticles.begin()->iam = store;
   }
}

/*
 * basic neighbours function. This time it loops over pInfos, which is sorted
 * spatially. As Thefunct is calculated over neighbouring particles pairs, this
 * will improve cache performance, as close particles have already been
 * accessed.
 */
template <void Thefunct(Cparticle &,Cparticle &,CglobalVars &)>
void CdataLL::neighbours() {
   for (vector<CpInfo>::iterator ppInfo = pInfos.begin();ppInfo!=pInfos.end();ppInfo++) {
      Cparticle *pp = ppInfo->p;
      if (ppInfo->neighbrs.empty()) {
         calcNeighbours(ppInfo->neighbrs,*pp);
         //cout << "Found "<<neighbrs[i].size()<<" neighbours for particle "<<i<<endl;
      }
      //for each neighbours run function
      for (vector<Cparticle *>::iterator pNeighbr = ppInfo->neighbrs.begin();pNeighbr!=ppInfo->neighbrs.end();pNeighbr++) {
         Thefunct(*pp,**pNeighbr,globals);
      }
   }
}

/*
 * same as above but adds an if function...
 */
template <void Thefunct(Cparticle &,Cparticle &,CglobalVars &),bool Iffunct(Cparticle &)>
void CdataLL::neighbours() {
   for (vector<CpInfo>::iterator ppInfo = pInfos.begin();ppInfo!=pInfos.end();ppInfo++) {
      Cparticle *pp = ppInfo->p;
      if (Iffunct(*pp)) {
         if (ppInfo->neighbrs.empty()) {
            //cout <<globals.mpiRank<<" before calc neigh"<<endl;
            calcNeighbours(ppInfo->neighbrs,*pp);
            //cout <<globals.mpiRank<<" after calc neigh"<<endl;
         }
         //for each neighbours run function
         for (vector<Cparticle *>::iterator pNeighbr = ppInfo->neighbrs.begin();pNeighbr!=ppInfo->neighbrs.end();pNeighbr++) {
            Thefunct(*pp,**pNeighbr,globals);
         }
      }
   }
}

template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &)>
void CdataLL::neighboursGroup() {
   for (vector<CpInfo>::iterator ppInfo = pInfos.begin();ppInfo!=pInfos.end();ppInfo++) {
      particleContainer::iterator pp = ppInfo->p;
      if (ppInfo->neighbrs.empty()) {
         calcNeighbours(ppInfo->neighbrs,*pp);
      }
      Thefunct(*pp,ppInfo->neighbrs,globals);
   }
}

template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &), bool Iffunct(Cparticle &)>
void CdataLL::neighboursGroup() {
   for (vector<CpInfo>::iterator ppInfo = pInfos.begin();ppInfo!=pInfos.end();ppInfo++) {
      Cparticle *pp = ppInfo->p;
      if (Iffunct(*pp)) {
         if (ppInfo->neighbrs.empty()) {
            calcNeighbours(ppInfo->neighbrs,*pp);
         }
         Thefunct(*pp,ppInfo->neighbrs,globals);
      }
   }
}


template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &), bool Iffunct(Cparticle &)>
void CdataLL::neighboursGroupAtRadius(const double radius) {
   const double radius2 = radius*radius;
   for (vector<CpInfo>::iterator ppInfo = pInfos.begin();ppInfo!=pInfos.end();ppInfo++) {
      Cparticle *pp = ppInfo->p;
      if (Iffunct(*pp)) {
         vector<Cparticle *> neighbrs;
         calcNeighboursAtRadius(neighbrs,*pp,radius,radius2);
         Thefunct(*pp,neighbrs,globals);
      }
   }
}

#ifdef LIQ_DEM
template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &), bool Iffunct(Cparticle &)>
void CdataLL::neighboursGroupCoupling() {
   for (vector<CpInfo>::iterator ppInfo = pInfos.begin();ppInfo!=pInfos.end();ppInfo++) {
      Cparticle *pp = ppInfo->p;
      if (Iffunct(*pp)) {
         vector<Cparticle *> neighbrs;
         calcNeighboursCoupling(neighbrs,*pp);
         Thefunct(*pp,neighbrs,globals);
      }
   }
}
#endif

template <class T, void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &,T &), bool Iffunct(Cparticle &)>
void CdataLL::neighboursGroup(T &param) {
   for (vector<CpInfo>::iterator ppInfo = pInfos.begin();ppInfo!=pInfos.end();ppInfo++) {
      Cparticle *pp = ppInfo->p;
      if (Iffunct(*pp)) {
         if (ppInfo->neighbrs.empty()) {
            calcNeighbours(ppInfo->neighbrs,*pp);
         }
         Thefunct(*pp,ppInfo->neighbrs,globals,param);
      }
   }
}

template <void Thefunct(Cparticle &,Cparticle &,CglobalVars &)>
void CdataLL::neighboursUsing(vector<Cparticle> &_ps) {
   vector<Cparticle *> _neighbrs;
   for (int i=0;i<_ps.size();i++) {
      _neighbrs.clear();
      calcNeighbours(_neighbrs,_ps[i]);
      //for each neighbours run function
      for (int j=0;j<_neighbrs.size();j++) {
         Cparticle *_pthep = _neighbrs[j];
         Thefunct(_ps[i],*_pthep,globals);
      }
   }
}

template <void Thefunct(Cparticle &,vector<Cparticle *> &,CglobalVars &)>
void CdataLL::neighboursUsing(vector<Cparticle> &_ps) {
   vector<Cparticle *> _neighbrs;
   int n = _ps.size();
   for (int i=0;i<n;i++) {
      _neighbrs.clear();
      calcNeighbours(_neighbrs,_ps[i]);
      Thefunct(_ps[i],_neighbrs,globals);
   }
}


template <void theFunct(Cparticle &,Cparticle &,CglobalVars &)>
void CdataLL::functionOverGrid(vector<Cparticle> &outPs,vect min,vect max,vectInt &gridDimN) {
   
   gridDimN = (max-min)/PSEP;
   vect newGridSep = 0.0;
   for (int i=0;i<NDIM;i++) {
      if (gridDimN[i]!=0) newGridSep[i] = (max[i]-min[i])/gridDimN[i];    
   }
   gridDimN += 1;
   int gridN = product(gridDimN); 
   outPs.resize(gridN);
  
   for (int i=0;i<gridDimN[0];i++) {
      if (NDIM > 1) {
         for (int j=0;j<gridDimN[1];j++) {
            int jOffset = j*gridDimN[0];
            if (NDIM > 2) {
               for (int k=0;k<gridDimN[2];k++) {
                  int kOffset = k*gridDimN[0]*gridDimN[1];
                  outPs[kOffset+jOffset+i].r = i*newGridSep[0],j*newGridSep[1],k*newGridSep[2];
                  outPs[kOffset+jOffset+i].r += min; 
               }
            } else {
               outPs[jOffset+i].r = i*newGridSep[0],j*newGridSep[1];
               outPs[jOffset+i].r += min; 
            }
         }
      } else {
         outPs[i].r = i*newGridSep[0];
         outPs[i].r += min; 
      } 
   }
   neighboursUsing<theFunct>(outPs);
}





template <class T, T readFunct(Cparticle &), void writeFunct(Cparticle &,T)>
void CdataLL::syncParticlesBetweenProcs() {
//NOTE: if the size sent per particle is greater than sizeof(CghostData) the send and recv buffers are too small!

   //fill up send buffers
   for (Array<int,NDIM>::iterator ap = globals.procNeighbrs.begin();ap != globals.procNeighbrs.end();ap++) {
      vectInt coords = ap.position();
      if (globals.procNeighbrs(coords)>=0) {
         T *p = (T *)syncBuffersSend(coords);
         for (int j=0;j<sendSizesGhosts(coords);j++) {
            *p = readFunct(*(ghostedParticles(coords)[j]));
            p++;
         }
      }
   }
  
   //send and recv everything (non-blocking)
   int sendRecvArraySize = 2*int(pow(3.0,NDIM));
   MPI_Request requestSendRecv[sendRecvArraySize];
   MPI_Status statusSendRecv[sendRecvArraySize];
   if (globals.mpiRank%2==0) {
      int upper = int(pow(3.0,NDIM));
      for (int i=0;i<upper;i++) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            sendRecvDataSync(coords,i,requestSendRecv+i*2,sizeof(T)*sendSizesGhosts(coords),syncBuffersSend(coords),sizeof(T)*recvSizesGhosts(coords),syncBuffersRecv(coords));
         }
         else {
            for (int j=0;j<2;j++) {
               requestSendRecv[i*2+j] = MPI_REQUEST_NULL;
            }
         }
      }
   } else {
      for (int i=int(pow(3.0,NDIM))-1;i>=0;i--) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            sendRecvDataSync(coords,i,requestSendRecv+i*2,sizeof(T)*sendSizesGhosts(coords),syncBuffersSend(coords),sizeof(T)*recvSizesGhosts(coords),syncBuffersRecv(coords));
         }
         else {
            for (int j=0;j<2;j++) {
               requestSendRecv[i*2+j] = MPI_REQUEST_NULL;
            }
         }
      } 
   }
   MPI_Waitall(sendRecvArraySize, requestSendRecv, statusSendRecv);   

   //sync new data with ghost particles
   vector<Cparticle>::iterator gp = procGhostParticles.begin();
   for (Array<int,NDIM>::iterator ap = globals.procNeighbrs.begin();ap != globals.procNeighbrs.end();ap++) {
      vectInt coords = ap.position();
      if (globals.procNeighbrs(coords)>= 0) {
         T *buffPtrT = (T *)syncBuffersRecv(coords);
         for (int z=0;z<recvSizesGhosts(coords);z++) {
            writeFunct(*gp,*buffPtrT);  
            buffPtrT++;
            gp++;
         }
      }
   }
}


template <class T, T readFunct(Cparticle &), void writeFunct(Cparticle &,T)>
void CdataLL::reverseSyncParticlesBetweenProcs() {
//NOTE: if the size sent per particle is greater than sizeof(CghostData) the send and recv buffers are too small!

   //fill up send buffers
   vector<Cparticle>::iterator gp = procGhostParticles.begin();
   for (Array<int,NDIM>::iterator ap = globals.procNeighbrs.begin();ap != globals.procNeighbrs.end();ap++) {
      vectInt coords = ap.position();
      if (globals.procNeighbrs(coords)>=0) {
         T *p = (T *)syncBuffersRecv(coords);
         for (int j=0;j<recvSizesGhosts(coords);j++) {
            *p = readFunct(*gp);
            p++;
            gp++;
         }
      }
   }
  
   //send and recv everything (non-blocking)
   int sendRecvArraySize = 2*int(pow(3.0,NDIM));
   MPI_Request requestSendRecv[sendRecvArraySize];
   MPI_Status statusSendRecv[sendRecvArraySize];
   if (globals.mpiRank%2==0) {
      int upper = int(pow(3.0,NDIM));
      for (int i=0;i<upper;i++) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            sendRecvDataSync(coords,i,requestSendRecv+i*2,sizeof(T)*recvSizesGhosts(coords),syncBuffersRecv(coords),sizeof(T)*sendSizesGhosts(coords),syncBuffersSend(coords));
         }
         else {
            for (int j=0;j<2;j++) {
               requestSendRecv[i*2+j] = MPI_REQUEST_NULL;
            }
         }
      }
   } else {
      for (int i=int(pow(3.0,NDIM))-1;i>=0;i--) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            sendRecvDataSync(coords,i,requestSendRecv+i*2,sizeof(T)*recvSizesGhosts(coords),syncBuffersRecv(coords),sizeof(T)*sendSizesGhosts(coords),syncBuffersSend(coords));
         }
         else {
            for (int j=0;j<2;j++) {
               requestSendRecv[i*2+j] = MPI_REQUEST_NULL;
            }
         }
      } 
   }
   MPI_Waitall(sendRecvArraySize, requestSendRecv, statusSendRecv);   

   //sync new data with ghost particles
   for (Array<int,NDIM>::iterator ap = globals.procNeighbrs.begin();ap != globals.procNeighbrs.end();ap++) {
      vectInt coords = ap.position();
      if (globals.procNeighbrs(coords)>= 0) {
         T *buffPtrT = (T *)syncBuffersSend(coords);
         for (int j=0;j<sendSizesGhosts(coords);j++) {
            writeFunct(*(ghostedParticles(coords)[j]),*buffPtrT);
            buffPtrT++;
         }
      }
   }
}
#endif
 
