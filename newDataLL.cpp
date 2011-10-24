#include "dataLL.impl.h"
#include "sort.h"

const unsigned int KEY_SIZE = 30;  // must be divisable by NDIM
const unsigned int KEY_DIM_BITS = KEY_SIZE/NDIM;
const unsigned int KEY_DIM_MAX = (unsigned int)floor(pow(2.0,(int)KEY_DIM_BITS)-1);

template <class particleT>
void CdataLL<particleT>::insertionSort(std::list<CpInfo> &v) {
   list<CpInfo>::iterator second = v.begin();
   second++;
   for (list<CpInfo>::iterator i = second;i!=v.end();i++) {
      for (list<CpInfo>::iterator j=v.begin();j!=i;j++) {
         if (j->key > i->key) {
            v.splice(i,v,j);
         }
      }
   }
}

template <class particleT>
void CdataLL<particleT>::sortPInfos() {
   for (vector<CpInfo>::iterator i=pInfos.begin()+1;i!=pInfos.end();i++) {
      for (vector<CpInfo>::iterator j=pInfos.begin();j!=i;j++) {
         if (j->key > i->key) {
            CpInfo tmpInfo = *j;
            *j = *i;
            j->ppInfoLink->ppInfo = &(*j);
            vector<CpInfo>::iterator jplus1 = j+1;
            for (vector<CpInfo >::iterator k=i;k!=jplus1;k--) {
               *k = *(k-1);
               k->ppInfoLink->ppInfo = &(*k);
            }
            *jplus1 = tmpInfo;
            jplus1->ppInfoLink->ppInfo = &(*jplus1);
         }
      }
   }
}

template <class particleT>
void CdataLL<particleT>::reset() {
   reset(1.0);
}

template <class particleT>
void CdataLL<particleT>::reset(const double inScale) {
   
   tagSortedParticlePointers.clear();

   scale = inScale;
   timeval t1,t2;
   gettimeofday(&t1,NULL);
   //removeGhosts();
   //insertGhosts();
   //neighbrs.resize(n);
   //cout << "Processor "<<globals.mpiRank<<": starting updateDomain..."<<endl;
   int psCapacity = ps.capacity();
   updateDomain();
   if ((ps.size() != pInfos.size())||(ps.size() != pInfoLinks.size())) {
      if (globals.sphStep != 0) {
         cout <<"Error: ps has "<<ps.size()<<" particles, pInfos has "<<pInfos.size()<<" particles, pInfoLinks has "<<pInfoLinks.size()<<" partcles"<<endl;
         exit(-1);
      }
   }
   if (ps.size() > psCapacity) {
      cerr << "Error: dataLL:reset(): ps has grown beyond its capacity! Bad things can happen......I'm outta here!"<<endl;
      exit(-1);
   }
   //cout << "Processor "<<globals.mpiRank<<": finished updateDomain..."<<endl;
   n = ps.size();
   calcDomainLimits();
   
   //for (Array<CsmVertex,NDIM>::iterator i=smGrid.begin();i!=smGrid.end();i++) {
   //   i->v = 0;
   //}
   //cout << "Processor "<<globals.mpiRank<<": rmax = "<<rmax<<" rmin = "<<rmin<<" ..."<<endl;
   cells.resize(static_cast<vectInt>((rmax-rmin)/(KERNAL_RADIUS*scale*hmax))+6);
   for (Array<vector<Cparticle *>,NDIM>::iterator i=cells.begin();i!=cells.end();i++) {
      i->clear();
   }
   //cout << "Processor "<<globals.mpiRank<<": ps has "<<ps.size()<<" particles. pInfos has "<<pInfos.size()<<" particle infos..."<<endl;
   for (vector<particleT>::iterator p = ps.begin();p != ps.end();p++) {
      //cout << "Processor "<<globals.mpiRank<<": inserting particle at r = "<<ppInfo->p->r<<endl;
      insert(*p,scale);
   }
   //cout << "Processor "<<globals.mpiRank<<": inserting "<<procGhostParticles.size()<<"ghost particles..."<<endl;
   for (vector<Cparticle>::iterator p=procGhostParticles.begin();p!= procGhostParticles.end();p++) {
      //cout << "Processor "<<globals.mpiRank<<": inserting ghost particle at r = "<<p->r<<endl;
      insert(*p,scale);
   }
   //cout << "Processor "<<globals.mpiRank<<": finsihed inserting ghost particles..."<<endl;
   gettimeofday(&t2,NULL);
   globals.wtTotalReset.tv_sec += t2.tv_sec-t1.tv_sec;
   globals.wtTotalReset.tv_usec += t2.tv_usec-t1.tv_usec;
}

template <class particleT>
void CdataLL<particleT>::insertNewParticle(Cparticle &newP) {
   ps.push_back(newP);
   n++;
   Cparticle *p = &(ps.back());
   CpInfo newpInfo;
   newpInfo.p = p;
   newpInfo.key = calcKey(newpInfo.p->r);
   pInfos.push_back(newpInfo);
   CpInfoLink newpInfoLink;
   newpInfoLink.ppInfo = &(pInfos.back());
   pInfoLinks.push_back(newpInfoLink);
   pInfos.back().ppInfoLink = &(pInfoLinks.back());
   insert(*p,1);
}

template <class particleT>
void CdataLL<particleT>::markForDeletion(Cparticle &p) {
   p.tag = -111;
}

template <class particleT>
void CdataLL::deleteParticles(particleT &p) {
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      particleT *thisP = ppInfo->p;
      if (thisP->tag == -111) {
         deleteParticle(ppInfoLink,ppInfo);
         continue;
      }
      ppInfoLink++;
   }
}

template <class particleT>
void CdataLL<particleT>::deleteParticle(vector<CpInfoLink>::iterator &ppInfoLink,CpInfo *ppInfo) {
   particleT *thisP = ppInfo->p;
   if (ppInfoLink != pInfoLinks.end()-1) {
      *thisP = ps.back();
      *ppInfoLink = pInfoLinks.back();

      ppInfoLink->ppInfo->p = thisP;
      ppInfoLink->ppInfo->ppInfoLink = &(*ppInfoLink);

      ps.pop_back();
      pInfoLinks.pop_back();
   } else {
      ps.pop_back();
      pInfoLinks.pop_back();
      ppInfoLink = pInfoLinks.end();
   }

   if (ppInfo != &(pInfos.back())) {
      *ppInfo = pInfos.back();
      ppInfo->ppInfoLink->ppInfo = ppInfo;
      pInfos.pop_back();
   } else {
      pInfos.pop_back();
   }
}



template <class particleT>
void CdataLL<particleT>::updateDomain() {
   //for each neighbouring proc, send real then ghost particles
   
   timeval t11,t22;
   gettimeofday(&t11,NULL);

   Array<int,NDIM> ghostIndicies(3);
   Array<int,NDIM> particleIndicies(3);
   ghostIndicies = 0;
   particleIndicies = 0;

   procGhostParticles.clear();

   //cout << "Processor "<<globals.mpiRank<<": creating dmin,dspace,dgmin,dgspace..."<<endl;
   vect dmin,dspace,dgmin,dgspace;
   for (int i=0;i<NDIM;i++) {
      dmin[i] = globals.procDomain[i*2];
      dspace[i] = globals.procDomain[i*2+1] - dmin[i];
      vectInt coords = 1;
      coords[i] = 0;
      dgmin[i] = dmin[i] + procGhostMax2H(coords);
      dgspace[i] = dspace[i] - procGhostMax2H(coords);
      //cout <<"procGhostMax2h = "<<procGhostMax2H(coords)<<endl;
      coords[i] = 2;
      dgspace[i] -= procGhostMax2H(coords);
      //cout <<"procGhostMax2h = "<<procGhostMax2H(coords)<<endl;
   }

   //cout << "Processor "<<globals.mpiRank<<": dmin = "<<dmin<<" dspace = "<<dspace<<" dgmin = "<<dgmin<<" dgspace = "<<dgspace<<endl;
   //cout << "Processor "<<globals.mpiRank<<": filling up buffers..."<<endl;
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      ppInfo->neighbrs.clear();
      ppInfo->scaledNeighbrs.clear();
      Cparticle *thisP = ppInfo->p;
      if (thisP->tag == -111) {
         deleteParticle(ppInfoLink,ppInfo);
         continue;
      }
      vectInt outCoords = static_cast<vectInt>((thisP->r-dmin)/dspace+1);
      if (globals.procNeighbrs(outCoords) >= 0) {
         //send particle to neighbour
         //add to real particle to send
         //cout << "Processor "<<globals.mpiRank<<": sending particle at "<<thisP->r<<" to processor "<<globals.procNeighbrs(outCoords)<<" ..."<<endl;
         if (particleIndicies(outCoords)>=pBufferSize) {
            cerr << "particleBuffers are full, exiting...."<<endl;
            exit(-1);
         }
         particleBuffersSend(outCoords)[particleIndicies(outCoords)] = *thisP;
         particleIndicies(outCoords)++;
         //delete from this particle array add to ghost array
         //ppInfo = pInfos.erase(ppInfo);
         //if (ppInfo != pInfos.begin()) ppInfo--;
         //procGhostParticles.push_back(*thisP);
         //procGhostParticles.back().iam = ghost;
         //ps.erase(thisP);
         deleteParticle(ppInfoLink,ppInfo);
         continue;
      }
      ppInfoLink++;
   }
   //cout << "Processor "<<globals.mpiRank<<": sending particle stuff..."<<endl;
   Array<int,NDIM> particleSizesSend(3);
   Array<int,NDIM> particleSizesRecv(3);
   particleSizesSend = particleIndicies*sizeof(Cparticle);
   int sendRecvArraySize = 4*int(pow(3.0,NDIM));
   MPI_Request requestSendRecv[sendRecvArraySize];
   MPI_Status statusSendRecv[sendRecvArraySize];

   timeval t1,t2;
   gettimeofday(&t1,NULL);
   if (globals.mpiRank%2==0) {
      int upper = int(pow(3.0,NDIM));
      for (int i=0;i<upper;i++) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            sendRecvParticles(coords,i,requestSendRecv+i*4,&(particleSizesSend(coords)),&(particleSizesRecv(coords))); 
         }
         else {
            for (int j=0;j<4;j++) {
               requestSendRecv[i*4+j] = MPI_REQUEST_NULL;
            }
         }
      }
   } else {
      for (int i=int(pow(3.0,NDIM))-1;i>=0;i--) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            sendRecvParticles(coords,i,requestSendRecv+i*4,&(particleSizesSend(coords)),&(particleSizesRecv(coords))); 
         }
         else {
            for (int j=0;j<4;j++) {
               requestSendRecv[i*4+j] = MPI_REQUEST_NULL;
            }
         }
      } 
   }
   MPI_Waitall(sendRecvArraySize, requestSendRecv, statusSendRecv);   

   gettimeofday(&t2,NULL);
   globals.wtTotalMPI.tv_sec += t2.tv_sec-t1.tv_sec;
   globals.wtTotalMPI.tv_usec += t2.tv_usec-t1.tv_usec;
     
 
   particleIndicies = particleSizesRecv/sizeof(Cparticle);

   //cout << "Processor "<<globals.mpiRank<<": finished sending bytes"<<endl;
   for (Array<int,NDIM>::iterator ap = globals.procNeighbrs.begin();ap != globals.procNeighbrs.end();ap++) {
      vectInt coords = ap.position();
      int neighbr = globals.procNeighbrs(coords);
      if (neighbr >= 0) {
         //put real particles in ps list
         double maxh = hmax;

         for (int z=0;z<particleIndicies(coords);z++) {
            Cparticle *recvP = &(particleBuffersRecv(coords)[z]);
            //cout << "Processor "<<globals.mpiRank<<": got a real particle at "<<recvP->r<<" ..."<<endl;
            ps.push_back(*recvP);
            Cparticle *p = &(ps.back());
            for (int i=0;i<NDIM;i++) {
               //if ((coords[i]==0)&&(p->r[i] > globals.procDomain[i*2])) {
               if ((PERIODIC[i])&&(coords[i]==0)&&(RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)) {
                  //cout << "Processor "<<globals.mpiRank<<": k = 0,  mapped a particle at r="<<p->r<<" to ";
                  p->r[i] = RMIN[i]-RMAX[i]+p->r[i]; 
                  //cout <<" about to drop dens by "<<DENS_DROP[i]<<endl;
                  if (p->iam==sph)
                     p->dens += DENS_DROP[i];
                  //cout <<" done "<<endl;
                  #ifdef SLK
                  p->currR[i] = RMIN[i]-RMAX[i]+p->currR[i];
                  #endif
                  //cout << "r = "<<p->r<<endl;
               }
               //if ((coords[i]==2)&&(p->r[i] < globals.procDomain[i*2+1])) {
               if ((PERIODIC[i])&&(coords[i]==2)&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP)) {
                  //cout << "Processor "<<globals.mpiRank<<": k = 1, mapped a particle at r="<<p->r<<" to ";
                  p->r[i] = RMAX[i]-RMIN[i]+p->r[i]; 
                  //cout <<" about to drop dens by "<<DENS_DROP[i]<<endl;
                  if (p->iam==sph)
                     p->dens -= DENS_DROP[i];
                  //cout <<" done "<<endl;
                  #ifdef SLK
                  p->currR[i] = RMAX[i]-RMIN[i]+p->currR[i];
                  #endif
                  //cout << "r = "<<p->r<<endl;
               }
               if ((GHOST[2*i])&&(coords[i]==0)&&(RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)) {
                  //cout << "Processor "<<globals.mpiRank<<": k = 0,  ghosted a particle at r="<<p.r<<" to ";
                  p->r[i] = 2.0*RMIN[i]-p->r[i];
                  if (GHOST[2*i]==1) {
                     p->v[i] = -p->v[i];
                     p->vhat[i] = -p->vhat[i];
                  } else {
                     p->v = -p->v;
                     p->vhat = -p->vhat;
                  }
                  //recurseGhost(p,i,procGhostParticles,hmax);
                  //cout << "r = "<<p.r<<endl;
               }
               //if ((coords[i]==2)&&(p.r[i] < globals.procDomain[i*2+1])) {
               if ((GHOST[2*i+1])&&(coords[i]==2)&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP)) {
                  //cout << "Processor "<<globals.mpiRank<<": k = 0,  ghosted a particle at r="<<p.r<<" to ";
                  p->r[i] = 2.0*RMAX[i]-p->r[i];
                  if (GHOST[2*i+1]==1) {
                     p->v[i] = -p->v[i];
                     p->vhat[i] = -p->vhat[i];
                  } else {
                     p->v = -p->v;
                     p->vhat = -p->vhat;
                  }
                  //recurseGhost(p,i,procGhostParticles,hmax);
                  //cout << "r = "<<p.r<<endl;
               }
            }
            CpInfo newpInfo;
            newpInfo.p = p;
            newpInfo.key = calcKey(newpInfo.p->r);
            pInfos.push_back(newpInfo);
            CpInfoLink newpInfoLink;
            newpInfoLink.ppInfo = &(pInfos.back());
            pInfoLinks.push_back(newpInfoLink);
            pInfos.back().ppInfoLink = &(pInfoLinks.back());
         }
      }
   }

   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();ppInfoLink++) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      ppInfo->neighbrs.clear();
      ppInfo->scaledNeighbrs.clear();
      Cparticle *thisP = ppInfo->p;
      vectInt inCoords = static_cast<vectInt>((thisP->r-dgmin)/dgspace+1);
      //cout << "Processor "<<globals.mpiRank<<": r = "<<thisP->r<<" outCoords = "<<outCoords<<" inCoords = "<<inCoords<<" neighbr = "<<endl;
      if (any(inCoords!=1)) {
         vectInt lowerBounds;
         vectInt upperBounds;
         for (int i=0;i<NDIM;i++) {
            if (inCoords[i]<=1) {
               upperBounds[i] = 1;
               lowerBounds[i] = inCoords[i];
            } else {
               upperBounds[i] = inCoords[i];
               lowerBounds[i] = 1;
            }
         }
         vectInt extent = upperBounds-lowerBounds+1;
         Array<int,NDIM> dummy(lowerBounds,extent);

         for (Array<int,NDIM>::iterator pD = dummy.begin();pD !=dummy.end();pD++) {
            vectInt dCoords = pD.position();
            if (globals.procNeighbrs(dCoords) >= 0) {
               //add to ghost particle array
               //cout << "Processor "<<globals.mpiRank<<": sending particle as ghost at "<<thisP->r<<" to processor "<<globals.procNeighbrs(dCoords)<<" ..."<<endl;
               if (ghostIndicies(dCoords)>=gBufferSize) {
                  cerr << "ghostBuffers are full, exiting...."<<endl;
                  exit(-1);
               }
               ghostBuffersSend(dCoords)[ghostIndicies(dCoords)] = *thisP;
               ghostedParticles(dCoords)[ghostIndicies(dCoords)] = thisP;
               ghostIndicies(dCoords)++;
            }
         }
      }
   }
   sendSizesGhosts = ghostIndicies;

   //cout << "Processor "<<globals.mpiRank<<": sending ghost stuff..."<<endl;
   Array<int,NDIM> ghostSizesSend(3);
   Array<int,NDIM> ghostSizesRecv(3);
   ghostSizesSend = ghostIndicies*sizeof(CghostData);
   //int sendRecvArraySize = 8*int(pow(3.0,NDIM));
   //MPI_Request requestSendRecv[sendRecvArraySize];
   //MPI_Status statusSendRecv[sendRecvArraySize];

   gettimeofday(&t1,NULL);
   if (globals.mpiRank%2==0) {
      int upper = int(pow(3.0,NDIM));
      for (int i=0;i<upper;i++) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            sendRecvGhosts(coords,i,requestSendRecv+i*4,&(ghostSizesSend(coords)),&(ghostSizesRecv(coords))); 
         }
         else {
            for (int j=0;j<4;j++) {
               requestSendRecv[i*4+j] = MPI_REQUEST_NULL;
            }
         }
      }
   } else {
      for (int i=int(pow(3.0,NDIM))-1;i>=0;i--) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            sendRecvGhosts(coords,i,requestSendRecv+i*4,&(ghostSizesSend(coords)),&(ghostSizesRecv(coords))); 
         }
         else {
            for (int j=0;j<4;j++) {
               requestSendRecv[i*4+j] = MPI_REQUEST_NULL;
            }
         }
      } 
   }
   MPI_Waitall(sendRecvArraySize, requestSendRecv, statusSendRecv);   

   gettimeofday(&t2,NULL);
   globals.wtTotalMPI.tv_sec += t2.tv_sec-t1.tv_sec;
   globals.wtTotalMPI.tv_usec += t2.tv_usec-t1.tv_usec;
   //cout << "statuses are: ";
   //for (int i=0;i<sendRecvArraySize;i++) {
   //   cout <<statusSendRecv[i].MPI_ERROR<<' ';
   //}
   //cout << endl;
      
   ghostIndicies = ghostSizesRecv/sizeof(CghostData);
   recvSizesGhosts = ghostIndicies;

   //int indexGhostParticles = procGhostParticles.size();

   //cout << "Processor "<<globals.mpiRank<<": finished sending bytes"<<endl;
   for (Array<int,NDIM>::iterator ap = globals.procNeighbrs.begin();ap != globals.procNeighbrs.end();ap++) {
      vectInt coords = ap.position();
      int neighbr = globals.procNeighbrs(coords);
      if (neighbr >= 0) {
         //put ghost particles in ghost vectors
         //work out ghost max h's
         double maxh = hmax;

         for (int z=0;z<ghostIndicies(coords);z++) {
            //cout << "Processor "<<globals.mpiRank<<": got a ghost particle at "<<ghostBuffersRecv(coords)[z].r<<" ..."<<endl;
            Cparticle p;
            p = ghostBuffersRecv(coords)[z];
            for (int i=0;i<NDIM;i++) {
               //if ((coords[i]==0)&&(p.r[i] > globals.procDomain[i*2])) {
               if ((PERIODIC[i])&&(coords[i]==0)&&(RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)) {
                  //cout << "Processor "<<globals.mpiRank<<": k = 0,  mapped a particle at r="<<p.r[i]<<" to ";
                  p.r[i] = RMIN[i]-RMAX[i]+p.r[i]; 
                  //cout <<" about to drop dens by "<<DENS_DROP[i]<<endl;
                  if (p.iam==sph)
                     p.dens += DENS_DROP[i];
                  //cout <<" done. PRB =  "<<PRB<<" i = "<<i<<endl;
                  #ifdef SLK
                  p.currR[i] = RMIN[i]-RMAX[i]+p.currR[i];
                  #endif
                  //cout << "r = "<<p.r[i]<<endl;
               }
               //if ((coords[i]==2)&&(p.r[i] < globals.procDomain[i*2+1])) {
               if ((PERIODIC[i])&&(coords[i]==2)&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP)) {
                  //cout << "Processor "<<globals.mpiRank<<": k = 1, mapped a particle at r="<<p.r[i]<<" to ";
                  p.r[i] = RMAX[i]-RMIN[i]+p.r[i]; 
                  if (p.iam==sph)
                     p.dens -= DENS_DROP[i];
                  #ifdef SLK
                  p.currR[i] = RMAX[i]-RMIN[i]+p.currR[i];
                  #endif
                  //cout << "r = "<<p.r[i]<<endl;
               }
               if ((GHOST[2*i])&&(coords[i]==0)&&(RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)) {
                  //cout << "Processor "<<globals.mpiRank<<": k = 0,  ghosted a particle at r="<<p.r<<" to ";
                  p.r[i] = 2.0*RMIN[i]-p.r[i];
                  if (GHOST[2*i]==1) {
                     p.v[i] = -p.v[i];
                     p.vhat[i] = -p.vhat[i];
                  } else {
                     p.v = -p.v;
                     p.vhat = -p.vhat;
                  }
               }
               //if ((coords[i]==2)&&(p.r[i] < globals.procDomain[i*2+1])) {
               if ((GHOST[2*i+1])&&(coords[i]==2)&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP)) {
                  //cout << "Processor "<<globals.mpiRank<<": k = 0,  ghosted a particle at r="<<p.r<<" to ";
                  p.r[i] = 2.0*RMAX[i]-p.r[i];
                  if (GHOST[2*i+1]==1) {
                     p.v[i] = -p.v[i];
                     p.vhat[i] = -p.vhat[i];
                  } else {
                     p.v = -p.v;
                     p.vhat = -p.vhat;
                  }
               }
            }
            if (p.h>maxh) maxh = p.h;
            procGhostParticles.push_back(p);
         }
         procGhostMax2H(coords) = KERNAL_RADIUS*maxh;
      }
   }

   //if (indexGhostParticles < procGhostParticles.size()) {
   //   startOfRecvGhostParticles = procGhostParticles.begin()+indexGhostParticles;
   //} else {
   //   startOfRecvGhostParticles = procGhostParticles.end();
   //}

   gettimeofday(&t22,NULL);
   globals.wtTotalUpdateDomain.tv_sec += t22.tv_sec-t11.tv_sec;
   globals.wtTotalUpdateDomain.tv_usec += t22.tv_usec-t11.tv_usec;
}

template <class particleT>
void CdataLL<particleT>::sendRecvData(const vectInt coords,const int num,MPI_Request *request,int *ghostSizesSend,int *ghostSizesRecv,int *particleSizesSend,int *particleSizesRecv) {
   int neighbr = globals.procNeighbrs(coords);
   vectInt middle = 1;
   vectInt split = 3;
   vectInt oppCoords = 2-coords;
   for (int i=0;i<NDIM;i++) {
      if (((coords[i]==2)&&(GHOST[2*i+1])&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP))||
          ((coords[i]==0)&&(GHOST[2*i])&&((RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)))) { 
         oppCoords[i] = coords[i];
      }
   }
   int oppNum = Nmisc::coordsToNum(oppCoords,split);
        
   MPI_Status status;
   //send ghost array size
   //cout << "Processor "<<globals.mpiRank<<": sending stuff to neighbour "<<neighbr<<" at coords: "<<coords<<endl;
   //cout << " neighbrs are "<<globals.procNeighbrs<<endl;
   //cout << "Processor "<<globals.mpiRank<<": sending "<<*ghostSizesSend<<" bytes"<<endl;
   MPI_Isend(ghostSizesSend,1,MPI_INT,neighbr,100+num,MPI_COMM_WORLD,request);
   MPI_Irecv(ghostSizesRecv,1,MPI_INT,neighbr,100+oppNum,MPI_COMM_WORLD,request+1);
   //send ghost array
   MPI_Isend(ghostBuffersSend(coords),*ghostSizesSend,MPI_BYTE,neighbr,200+num,MPI_COMM_WORLD,request+2);
   MPI_Irecv(ghostBuffersRecv(coords),gBufferSize*sizeof(CghostData),MPI_BYTE,neighbr,200+oppNum,MPI_COMM_WORLD,request+3);

   //send real array 
   //cout << "Processor "<<globals.mpiRank<<": sending "<<*particleSizesSend<<" bytes"<<endl;
   MPI_Isend(particleSizesSend,1,MPI_INT,neighbr,300+num,MPI_COMM_WORLD,request+4);
   MPI_Irecv(particleSizesRecv,1,MPI_INT,neighbr,300+oppNum,MPI_COMM_WORLD,request+5);
   MPI_Isend(particleBuffersSend(coords),*particleSizesSend,MPI_BYTE,neighbr,400+num,MPI_COMM_WORLD,request+6);
   MPI_Irecv(particleBuffersRecv(coords),pBufferSize*sizeof(Cparticle),MPI_BYTE,neighbr,400+oppNum,MPI_COMM_WORLD,request+7);
}

template <class particleT>
void CdataLL<particleT>::sendRecvParticles(const vectInt coords,const int num,MPI_Request *request,int *particleSizesSend,int *particleSizesRecv) {
   int neighbr = globals.procNeighbrs(coords);
   vectInt middle = 1;
   vectInt split = 3;
   vectInt oppCoords = 2-coords;
   for (int i=0;i<NDIM;i++) {
      if (((coords[i]==2)&&(GHOST[2*i+1])&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP))||
          ((coords[i]==0)&&(GHOST[2*i])&&((RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)))) { 
         oppCoords[i] = coords[i];
      }
   }
   int oppNum = Nmisc::coordsToNum(oppCoords,split);
        
   MPI_Status status;
   //cout << "Processor "<<globals.mpiRank<<": sending stuff to neighbour "<<neighbr<<" at coords: "<<coords<<endl;
   //cout << " neighbrs are "<<globals.procNeighbrs<<endl;
   //cout << "Processor "<<globals.mpiRank<<": sending "<<*ghostSizesSend<<" bytes"<<endl;

   //send real array 
   //cout << "Processor "<<globals.mpiRank<<": sending "<<*particleSizesSend<<" bytes"<<endl;
   MPI_Isend(particleSizesSend,1,MPI_INT,neighbr,300+num,MPI_COMM_WORLD,request);
   MPI_Irecv(particleSizesRecv,1,MPI_INT,neighbr,300+oppNum,MPI_COMM_WORLD,request+1);
   MPI_Isend(particleBuffersSend(coords),*particleSizesSend,MPI_BYTE,neighbr,400+num,MPI_COMM_WORLD,request+2);
   MPI_Irecv(particleBuffersRecv(coords),pBufferSize*sizeof(Cparticle),MPI_BYTE,neighbr,400+oppNum,MPI_COMM_WORLD,request+3);
}

template <class particleT>
void CdataLL<particleT>::sendRecvGhosts(const vectInt coords,const int num,MPI_Request *request,int *ghostSizesSend,int *ghostSizesRecv) {
   int neighbr = globals.procNeighbrs(coords);
   vectInt middle = 1;
   vectInt split = 3;
   vectInt oppCoords = 2-coords;
   for (int i=0;i<NDIM;i++) {
      if (((coords[i]==2)&&(GHOST[2*i+1])&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP))||
          ((coords[i]==0)&&(GHOST[2*i])&&((RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)))) { 
         oppCoords[i] = coords[i];
      }
   }
   int oppNum = Nmisc::coordsToNum(oppCoords,split);
        
   MPI_Status status;
   //send ghost array size
   //cout << "Processor "<<globals.mpiRank<<": sending stuff to neighbour "<<neighbr<<" at coords: "<<coords<<endl;
   //cout << " neighbrs are "<<globals.procNeighbrs<<endl;
   //cout << "Processor "<<globals.mpiRank<<": sending "<<*ghostSizesSend<<" bytes"<<endl;
   MPI_Isend(ghostSizesSend,1,MPI_INT,neighbr,100+num,MPI_COMM_WORLD,request);
   MPI_Irecv(ghostSizesRecv,1,MPI_INT,neighbr,100+oppNum,MPI_COMM_WORLD,request+1);
   //send ghost array
   MPI_Isend(ghostBuffersSend(coords),*ghostSizesSend,MPI_BYTE,neighbr,200+num,MPI_COMM_WORLD,request+2);
   MPI_Irecv(ghostBuffersRecv(coords),gBufferSize*sizeof(CghostData),MPI_BYTE,neighbr,200+oppNum,MPI_COMM_WORLD,request+3);
}

template <class particleT>
void CdataLL<particleT>::sendRecvDataSync(const vectInt coords,const int num,MPI_Request *request,int sendSize, void *sendBuffer, int recvSize, void *recvBuffer) {
   int neighbr = globals.procNeighbrs(coords);
   vectInt middle = 1;
   vectInt split = 3;
   vectInt oppCoords = 2-coords;
   for (int i=0;i<NDIM;i++) {
      if (((coords[i]==2)&&(GHOST[2*i+1])&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP))||
          ((coords[i]==0)&&(GHOST[2*i])&&((RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)))) { 
         oppCoords[i] = coords[i];
      }
   }
   int oppNum = Nmisc::coordsToNum(oppCoords,split);
        
   MPI_Status status;
   //send ghost array size
   //cout << "Processor "<<globals.mpiRank<<": sending stuff to neighbour "<<neighbr<<" at coords: "<<coords<<endl;
   //cout << " neighbrs are "<<globals.procNeighbrs<<endl;
   //cout << "Processor "<<globals.mpiRank<<": sending "<<*ghostSizesSend<<" bytes"<<endl;
   //send ghost array
   MPI_Isend(sendBuffer,sendSize,MPI_BYTE,neighbr,500+num,MPI_COMM_WORLD,request);
   MPI_Irecv(recvBuffer,recvSize,MPI_BYTE,neighbr,500+oppNum,MPI_COMM_WORLD,request+1);
}

template <class particleT>
void CdataLL<particleT>::setGlobalTimestep(double dt) {
   double sendbuf = dt;
   double recvbuf[globals.mpiSize];
   //cout << "Processor "<<globals.mpiRank<<": starting setGlobalTimestep...sending dt = "<<dt<<endl;
   //for (int i=0;i<globals.mpiSize;i++) {
   //   sendbuf[i] = dt;
   //}
   MPI_Allgather(&sendbuf, 1, MPI_DOUBLE, recvbuf, 1, MPI_DOUBLE, MPI_COMM_WORLD); 
   globals.dt = dt;
   for (int i=0;i<globals.mpiSize;i++) {
      if (recvbuf[i]<globals.dt) globals.dt = recvbuf[i];
   }
   //cout << "Processor "<<globals.mpiRank<<": finished setGlobalTimestep...dt = "<<globals.dt<<endl;
}


template <class particleT>
void CdataLL<particleT>::sumOverProcs(vect *data,int num) {
//should make a general one: function over procs
   
   int newNum = num*sizeof(vect);
   vect *sendbuf = data;
   vect recvbuf[globals.mpiSize*num];
   //cout << "Processor "<<globals.mpiRank<<": starting sumOverProcs...sending data = "<<endl;
   //for (int i=0;i<num;i++) {
   //   cout << data[i]<<" ";
   //}
   //cout <<endl;

   MPI_Gather(sendbuf,newNum,MPI_BYTE,recvbuf,newNum,MPI_BYTE,0,MPI_COMM_WORLD);
   if (globals.mpiRank==0) {
      for (int i=1;i<globals.mpiSize;i++) {
         for (int j=0;j<num;j++) {
            recvbuf[j] += recvbuf[i*num+j];
         }
      }
      for (int i=1;i<globals.mpiSize;i++) {
         for (int j=0;j<num;j++) {
            recvbuf[i*num+j] = recvbuf[j];
         }
      }
   }
   MPI_Scatter(recvbuf,newNum,MPI_BYTE,sendbuf,newNum,MPI_BYTE,0,MPI_COMM_WORLD);
   //cout << "Processor "<<globals.mpiRank<<": finishing sumOverProcs...recv data = "<<endl;
   //for (int i=0;i<num;i++) {
   //   cout << data[i]<<" ";
   //}
   //cout <<endl;
}

template <class particleT>
void CdataLL<particleT>::sumOverProcs(double *data,int num) {
//should make a general one: function over procs
   
   int newNum = num*sizeof(double);
   double *sendbuf = data;
   double recvbuf[globals.mpiSize*num];
   //cout << "Processor "<<globals.mpiRank<<": starting sumOverProcs...sending data = "<<endl;
   //for (int i=0;i<num;i++) {
   //   cout << data[i]<<" ";
   //}
   //cout <<endl;

   MPI_Gather(sendbuf,newNum,MPI_BYTE,recvbuf,newNum,MPI_BYTE,0,MPI_COMM_WORLD);
   if (globals.mpiRank==0) {
      for (int i=1;i<globals.mpiSize;i++) {
         for (int j=0;j<num;j++) {
            recvbuf[j] += recvbuf[i*num+j];
         }
      }
      for (int i=1;i<globals.mpiSize;i++) {
         for (int j=0;j<num;j++) {
            recvbuf[i*num+j] = recvbuf[j];
         }
      }
   }
   MPI_Scatter(recvbuf,newNum,MPI_BYTE,sendbuf,newNum,MPI_BYTE,0,MPI_COMM_WORLD);
   //cout << "Processor "<<globals.mpiRank<<": finishing sumOverProcs...recv data = "<<endl;
   //for (int i=0;i<num;i++) {
   //   cout << data[i]<<" ";
   //}
   //cout <<endl;
}

template <class particleT>
void CdataLL<particleT>::syncSmGrid() {
   //Array<CsmVertex*,NDIM> sendBufs(3);      
   //Array<CsmVertex*,NDIM> recvBufs(3);
   //Array<int,NDIM> sizes(3);
   //for (Array<CsmVertex*,NDIM>::iterator i = sendBufs.begin();i != sendBufs.end();i++) {
   for (Array<CsmVertex*,NDIM>::iterator i = smGridBuffersSend.begin();i != smGridBuffersSend.end();i++) {
      vectInt coords = i.position();
      int neighbr = globals.procNeighbrs(coords);
      if (neighbr >= 0) {
        /* vectInt upperBounds,lowerBounds;
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
         Array<CsmVertex,NDIM> smGridView = smGrid(subdomain);
         *i = new CsmVertex[smGridView.size()];
         recvBufs(coords) = new CsmVertex[smGridView.size()];
         sizes(coords) = smGridView.size();
      */

         int j=0;
         for (Array<CsmVertex,NDIM>::iterator k=smGridViews(coords).begin();k!=smGridViews(coords).end();k++) {
            (*i)[j] = *k;
            j++;
         }
      }
   }
          
   int sendRecvArraySize = 2*int(pow(3.0,NDIM));
   MPI_Request requestSendRecv[sendRecvArraySize];
   MPI_Status statusSendRecv[sendRecvArraySize];

   if (globals.mpiRank%2==0) {
      int upper = int(pow(3.0,NDIM));
      for (int i=0;i<upper;i++) {
         vectInt split = 3;
         vectInt coords = Nmisc::numToCoords(i,split);
         if (globals.procNeighbrs(coords)>=0) {
            //sendRecvSmGridData(coords,i,requestSendRecv+i*2,sendBufs(coords),recvBufs(coords),sizes(coords));
            sendRecvSmGridData(coords,i,requestSendRecv+i*2,smGridBuffersSend(coords),smGridBuffersRecv(coords),smGridViews(coords).size());
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
            //sendRecvSmGridData(coords,i,requestSendRecv+i*2,sendBufs(coords),recvBufs(coords),sizes(coords));
            sendRecvSmGridData(coords,i,requestSendRecv+i*2,smGridBuffersSend(coords),smGridBuffersRecv(coords),smGridViews(coords).size());
         }
         else {
            for (int j=0;j<2;j++) {
               requestSendRecv[i*2+j] = MPI_REQUEST_NULL;
            }
         }
      } 
   }
   MPI_Waitall(sendRecvArraySize, requestSendRecv, statusSendRecv);   

   //for (Array<CsmVertex*,NDIM>::iterator i = recvBufs.begin();i != recvBufs.end();i++) {
   for (Array<CsmVertex*,NDIM>::iterator i = smGridBuffersRecv.begin();i != smGridBuffersRecv.end();i++) {
      vectInt coords = i.position();
      int neighbr = globals.procNeighbrs(coords);
      if (neighbr >= 0) {
         /*vectInt upperBounds,lowerBounds;
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
         //cout <<"Processor "<<globals.mpiRank<<": coords = "<<coords<<" lowerbounds = "<<lowerBounds<<" upperbounds = "<<upperBounds<<endl;
         RectDomain<NDIM> subdomain(lowerBounds, upperBounds);
         Array<CsmVertex,NDIM> smGridView = smGrid(subdomain);
         */
         int j=0;
         for (Array<CsmVertex,NDIM>::iterator k=smGridViews(coords).begin();k!=smGridViews(coords).end();k++) {
            if (any(abs(k->pos-(*i)[j].pos)>=SMOOTH_GRID_SIZE/2.0)) {
               cerr << "Processor "<<globals.mpiRank<<" error: syncing wrong vertecies with each other! read "<<(*i)[j].pos<<" write "<<k->pos << endl;
               exit(-1);
            }
            k->v += (*i)[j].v; 
            //k->v = 1;
            k->num += (*i)[j].num; 
            //cout <<"Grid point at pos = "<<k->pos<<" num after = "<<k->num<<" just added "<<(*i)[j].num<<endl;
            j++;
         }
         //delete [] (*i);
         //delete [] sendBufs(coords);
      }
   }
}

template <class particleT>
void CdataLL<particleT>::sendRecvSmGridData(vectInt coords,int num,MPI_Request *request,CsmVertex *sendBuf,CsmVertex *recvBuf,int size) {
   int neighbr = globals.procNeighbrs(coords);
   vectInt middle = 1;
   vectInt split = 3;
   vectInt oppCoords = 2-coords;
   int oppNum = Nmisc::coordsToNum(oppCoords,split);
        
   MPI_Status status;
   //send and recv
   MPI_Isend(sendBuf,size*sizeof(CsmVertex),MPI_BYTE,neighbr,600+num,MPI_COMM_WORLD,request);
   MPI_Irecv(recvBuf,size*sizeof(CsmVertex),MPI_BYTE,neighbr,600+oppNum,MPI_COMM_WORLD,request+1);
}

template <class particleT>
void CdataLL<particleT>::getVerticesFromPosition(vect pos,vectInt &lowerB) {
   vect scalePos;
   for (int i=0;i<NDIM;i++) {
      scalePos[i] = (pos[i]-smGridMin[i])*(1.0/SMOOTH_GRID_SIZE);
      if (floor(scalePos[i])==ceil(scalePos[i])) {
         if (scalePos[i]!=0.0) {
            scalePos[i] -= SMOOTH_GRID_SIZE*0.001;
         } else {
            scalePos[i] += SMOOTH_GRID_SIZE*0.001;
         }
      }
   }
   //RectDomain<NDIM> subdomain(floor(scalePos), ceil(scalePos));
   //Array<CsmVertex,NDIM> test = smGrid(subdomain);
   //cout <<"Processor "<<globals.mpiRank<<": getting vertices. pos = "<<pos<<" scalePos = "<<scalePos<<" smGridMin = "<<smGridMin<<" minpos = "<<test(0,0).pos<<" maxpos = "<<test(1,1).pos<<endl;
   lowerB = floor(scalePos);
}

template <class particleT>
unsigned long CdataLL<particleT>::calcKey(vect r) {
   unsigned int rint[NDIM];
   unsigned long key = 0;
   //cout << "converting r = "<<r<<" into key..."<<endl;
   for (int i=0;i<NDIM;i++) {
      rint[i] = (unsigned int)floor((r[i]-rmin[i])*KEY_DIM_MAX/(rmax[i]-rmin[i]));
      //cout << "integer from dim "<<i<<" is "<<rint[i]<<" (rmax = "<<rmax<<" rmin = "<<rmin<<" KEY_DIM_MAX = "<<KEY_DIM_MAX<<" )"<<endl;
   }
   unsigned int mask = 0x01<<(KEY_DIM_BITS-1);
   //cout << "KEY_DIM_BITS = "<<KEY_DIM_BITS<<endl;
   for (int i=0;i<KEY_DIM_BITS;i++) {
      for (int j=0;j<NDIM;j++) {
         key <<= 1;
         key |= (rint[j] & mask) >> (KEY_DIM_BITS-1-i);
      }
      mask >>= 1;
   } 
   //cout << "calculating key from r="<<r<<". Result is key="<<key<<endl;
   return key; 
}

template <class particleT>
void CdataLL<particleT>::updateParticles() {
   initPinfos();
   checkPinfos();
}

template <class particleT>
void CdataLL<particleT>::initPinfos() {
   pInfos.clear();
   pInfoLinks.clear();
   pInfos.reserve(ps.capacity());
   pInfoLinks.reserve(ps.capacity());
   for (vector<particleT>::iterator thisP = ps.begin();thisP!=ps.end();thisP++) {
      CpInfo newpInfo;
      newpInfo.p = &(*thisP);
      newpInfo.key = calcKey(thisP->r);
      pInfos.push_back(newpInfo);

      CpInfoLink newpInfoLink;
      newpInfoLink.ppInfo = &(pInfos.back());
      pInfoLinks.push_back(newpInfoLink);
      pInfos.back().ppInfoLink = &(pInfoLinks.back());

      vector<CpInfo>::iterator ppInfo = pInfos.end()-1;
      vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.end()-1;
   }
}

template <class particleT>
void CdataLL<particleT>::checkPinfos() {
   vector<particleT>::iterator p = ps.begin();
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();ppInfoLink++) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      if (ppInfo->ppInfoLink != &(*ppInfoLink)) {
         cerr <<"Error: CdataLL<particleT>::initTraverseOrder: ppInfo and ppInfoLink are not consistant!"<<endl;
         exit(-1);
      }
      if (ppInfo->p != &(*p)) {
         cerr <<"Error: CdataLL<particleT>::initTraverseOrder: ppInfo, ppInfoLink and p are not consistant!"<<endl;
         exit(-1);
      }
      p++;
   }
}

template <class particleT>
void CdataLL<particleT>::initTraverseOrder() {
   //cout << "Processor "<<globals.mpiRank<<": starting initTraverseOrder..."<<endl;
   initPinfos();   
   sortPInfos();
   checkPinfos();
   //cout << "Processor "<<globals.mpiRank<<": finished initTraverseOrder..."<<endl;
}

template <class particleT>
void CdataLL<particleT>::calcTraverseOrder() {
   //cout << "Processor "<<globals.mpiRank<<": starting calcTraverseOrder..."<<endl;
   vector<particleT>::iterator p = ps.begin();
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();ppInfoLink++) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      if (ppInfo->ppInfoLink != &(*ppInfoLink)) {
         cerr <<"Error: CdataLL<particleT>::calcTraverseOrder: ppInfo and ppInfoLink are not consistant!"<<endl;
         exit(-1);
      }
      if (ppInfo->p != &(*p)) {
         cerr <<"Error: CdataLL<particleT>::calcTraverseOrder: ppInfo, ppInfoLink and p are not consistant!"<<endl;
         exit(-1);
      }
      ppInfo->key = calcKey(ppInfo->p->r);
      ppInfo->neighbrs.clear();
      p++;
   }
   sortPInfos();
   //cout << "Processor "<<globals.mpiRank<<": starting calcTraverseOrder..."<<endl;
}

template <class particleT>
void CdataLL<particleT>::calcDomainLimits() {
   for (int i=0;i<NDIM;i++) {
      rmin[i] = globals.procDomain[i*2];
      rmax[i] = globals.procDomain[i*2+1];
   }
   hmax = -100000;
   for (vector<particleT>::iterator thisP = ps.begin();thisP!=ps.end();thisP++) {
      for (int j=0;j<NDIM;j++) {
         if (thisP->r[j] < rmin[j]) { rmin[j] = thisP->r[j]; }
         if (thisP->r[j] > rmax[j]) { rmax[j] = thisP->r[j]; }
      }
     
      if (thisP->h > hmax) hmax = thisP->h;
   }
   //cout <<"found hmax = "<<hmax<<endl;
   for (int i=0;i<NDIM;i++) {
      vectInt coords = 1;
      coords[i] = 0;
      if (procGhostMax2H(coords)<hmax) procGhostMax2H(coords)=hmax;
      coords[i] = 2;
      if (procGhostMax2H(coords)<hmax) procGhostMax2H(coords)=hmax;
   }
   /*for (vector<Cparticle>::iterator thisP=procGhostParticles.begin();thisP!= procGhostParticles.end();thisP++) {
      for (int j=0;j<NDIM;j++) {
         if (thisP->r[j] < rmin[j]) { rmin[j] = thisP->r[j]; }
         if (thisP->r[j] > rmax[j]) { rmax[j] = thisP->r[j]; }
      }
     
      if (thisP->h > hmax) hmax = thisP->h;
   }*/
}


template <class particleT>
void CdataLL<particleT>::insertGhosts() {
   //TODO: this function is outdated, redo!
   /*Cparticle newP;
   for (list<Cparticle>::iterator thisP = ps.begin();thisP!=ps.end();thisP++) {
      for (int j=0;j<NDIM;j++) {
         if (thisP->r[j] < PMIN[j]) thisP->r[j] = PMAX[j]-PMIN[j]+thisP->r[j];
         if (thisP->r[j] > PMAX[j]) thisP->r[j] = PMIN[j]-PMAX[j]+thisP->r[j];
         if (thisP->r[j] < PMIN[j]+KERNAL_RADIUS*hmax) {
            newP = *thisP;
            newP.r[j] = PMAX[j]-PMIN[j]+thisP->r[j];
            if (newP.iam==sph) newP.iam = ghost;
            ps.push_back(newP);
         } 
         if (thisP->r[j] > PMAX[j]-KERNAL_RADIUS*hmax) {
            newP = *thisP;
            newP.r[j] = PMIN[j]-PMAX[j]+thisP->r[j];
            if (newP.iam==sph) newP.iam = ghost;
            ps.push_back(newP);
         }
      }
   }
   //ghostStart = ps.begin() + n;
   //ghostEnd = ps.end()-1;
   n = ps.size();
   */
}

template <class particleT>
void CdataLL<particleT>::removeGhosts() {
   //TODO: this function is outdated, redo!
   /*if (ghostStart!=ghostEnd) {
      ps.erase(ghostStart,ghostEnd+1);
   }
   n = ps.size();
   */
}
            

template <class particleT>
void CdataLL<particleT>::insert(Cparticle &p,const double scale) {
   vectInt cellI = getCellI(p.r);
   //cout << "Processor "<<globals.mpiRank<<": inserting particle in cell = "<<cellI<<". size of cell = "<<cells(cellI).size()<<endl;
   //for (int i=0;i<NDIM;i++) {
   //   if ((cellI[i]<0)||(cellI[i]>=cells.extent(i))) cout <<"error: cellI = "<<cellI<<endl;
   //} 

   cells(cellI).push_back(&p);
}

template <class particleT>
void CdataLL<particleT>::addIfNeighbr(Cparticle &_p,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd) {
   for (int i=0;i<_toAdd.size();i++) {
      if (_p.is_neighbr(*(_toAdd[i]))) {
         _neighbrs.push_back(_toAdd[i]);
      }
   }
}


template <class particleT>
void CdataLL<particleT>::calcNeighbours(vector<Cparticle *> &_neighbrs,Cparticle &_p) {
   vectInt cellI = getCellI(_p.r);
   vectInt tcellI = cellI;
   for (int i=cellI[0]-1;i<=cellI[0]+1;i++) {
      tcellI[0] = i;
      if (NDIM > 1) {
         for (int j=cellI[1]-1;j<=cellI[1]+1;j++) {
            tcellI[1] = j;
            if (NDIM > 2) {
               for (int k=cellI[2]-1;k<=cellI[2]+1;k++) {
                  tcellI[2] = k;
                  addIfNeighbr(_p,_neighbrs,cells(tcellI));
               }
            } else {
               addIfNeighbr(_p,_neighbrs,cells(tcellI));
            }
         }
      } else {
         addIfNeighbr(_p,_neighbrs,cells(tcellI));
      } 
   }
}


template <class particleT>
void CdataLL<particleT>::addIfBothNeighbr(Cparticle &_p,vector<Cparticle *> &_scaledNeighbrs,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd) {
   for (int i=0;i<_toAdd.size();i++) {
      if (_p.isScaledNeighbr(*(_toAdd[i]),scale)) {
         _scaledNeighbrs.push_back(_toAdd[i]);
         if (_p.is_neighbr(*(_toAdd[i]))) {
            _neighbrs.push_back(_toAdd[i]);
         }
      }
   }
}


template <class particleT>
void CdataLL<particleT>::calcBothNeighbours(vector<Cparticle *> &_scaledNeighbrs,vector<Cparticle *> &_neighbrs,Cparticle &_p) {
   vectInt cellI = getCellI(_p.r);
   vectInt tcellI = cellI;
   for (int i=cellI[0]-1;i<=cellI[0]+1;i++) {
      tcellI[0] = i;
      if (NDIM > 1) {
         for (int j=cellI[1]-1;j<=cellI[1]+1;j++) {
            tcellI[1] = j;
            if (NDIM > 2) {
               for (int k=cellI[2]-1;k<=cellI[2]+1;k++) {
                  tcellI[2] = k;
                  addIfBothNeighbr(_p,_scaledNeighbrs,_neighbrs,cells(tcellI));
               }
            } else {
               addIfBothNeighbr(_p,_scaledNeighbrs,_neighbrs,cells(tcellI));
            }
         }
      } else {
         addIfBothNeighbr(_p,_scaledNeighbrs,_neighbrs,cells(tcellI));
      } 
   }
}

template <class particleT>
void CdataLL<particleT>::addIfScaledNeighbr(Cparticle &_p,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd) {
   for (int i=0;i<_toAdd.size();i++) {
      if (_p.isScaledNeighbr(*(_toAdd[i]),scale)) {
         _neighbrs.push_back(_toAdd[i]);
      }
   }
}

template <class particleT>
void CdataLL<particleT>::calcScaledNeighbours(vector<Cparticle *> &_neighbrs,Cparticle &_p) {
   vectInt cellI = getCellI(_p.r);
   vectInt tcellI = cellI;
   for (int i=cellI[0]-1;i<=cellI[0]+1;i++) {
      tcellI[0] = i;
      if (NDIM > 1) {
         for (int j=cellI[1]-1;j<=cellI[1]+1;j++) {
            tcellI[1] = j;
            if (NDIM > 2) {
               for (int k=cellI[2]-1;k<=cellI[2]+1;k++) {
                  tcellI[2] = k;
                  addIfScaledNeighbr(_p,_neighbrs,cells(tcellI));
               }
            } else {
               addIfScaledNeighbr(_p,_neighbrs,cells(tcellI));
            }
         }
      } else {
         addIfScaledNeighbr(_p,_neighbrs,cells(tcellI));
      } 
   }
}


template <class particleT>
Cparticle *CdataLL<particleT>::getParticleFromTag(int tag) {
   //cout <<"get particle with tag = "<<tag<<endl;
   if (tagSortedParticlePointers.empty()) {
      cout <<"data struct empty, filling tagSortedParticlePointers"<<endl;
      tagSortedParticlePointers.resize(ps.size());
      for (vector<Cparticle>::reverse_iterator i = ps.rbegin();i!=ps.rend();i++) {
         if (tagSortedParticlePointers.size()<i->tag) {
            tagSortedParticlePointers.resize(i->tag*2);
         }
         tagSortedParticlePointers[i->tag-1].tag = i->tag;
         tagSortedParticlePointers[i->tag-1].p = &(*i);
      }
   }
   if (tagSortedParticlePointers[tag-1].tag != tag) {
      cerr << "Error in CdataLL<particleT>::getParticleFromTag(): the tag aint here!!!"<<endl;
      exit(-1);
   }
   return tagSortedParticlePointers[tag-1].p;
}

template <class particleT>
CpInfo *CdataLL<particleT>::getPinfoFromTag(int tag) {
   //cout <<"get particle with tag = "<<tag<<endl;
   if (tagSortedPinfoPointers.empty()) {
      cout <<"data struct empty, filling tagSortedPinfoPointers"<<endl;
      tagSortedPinfoPointers.resize(pInfos.size());
      for (vector<CpInfo>::reverse_iterator i = pInfos.rbegin();i!=pInfos.rend();i++) {
         if (tagSortedPinfoPointers.size()<i->p->tag) {
            tagSortedPinfoPointers.resize(i->p->tag*2);
         }
         tagSortedPinfoPointers[i->p->tag-1].tag = i->p->tag;
         tagSortedPinfoPointers[i->p->tag-1].pInfo = &(*i);
      }
   }
   if (tagSortedPinfoPointers[tag-1].tag != tag) {
      cerr << "Error in CdataLL<particleT>::getPinfoFromTag(): the tag aint here!!!"<<endl;
      exit(-1);
   }
   return tagSortedPinfoPointers[tag-1].pInfo;
}

int CdataLL<particleT>::calcBufferSize(vectInt coords) {
   int bufferSize = 1.0;
   for (int i=0;i<NDIM;i++) {
      if (coords[i]==0) {
         bufferSize *= int((RMAX[i]-RMIN[i])/PSEP);
      } else {
         bufferSize *= int(2.0*H/PSEP);
      }
   }
   return 2*bufferSize;
}    

    

         
   
   
