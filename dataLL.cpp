#include "dataLL.impl.h"
#include "sort.h"

const unsigned int KEY_SIZE = 30;  // must be divisable by NDIM
const unsigned int KEY_DIM_BITS = KEY_SIZE/NDIM;
const unsigned int KEY_DIM_MAX = (unsigned int)floor(pow(2.0,(int)KEY_DIM_BITS)-1);

CdataLL::CdataLL(particleContainer &in_ps, CglobalVars &in_globals,bool opt): ps(in_ps),globals(in_globals) {
   //initialise hmax, the maximum smoothing length
   hmax = H;
   
   cout << "size of particle array = "<<ps.capacity()*sizeof(Cparticle)/1024/1024<<" MB"<<endl;

   // allocate buffers used for mpi communication
   allocateBuffers();
         
   //keep track of hmax for each neighbouring cpu
   procGhostMax2H.resize(3);
   procGhostMax2H = 1.1*KERNAL_RADIUS*H;

   cout << "Processor "<<globals.mpiRank<<" of "<<globals.mpiSize<<" has domain ";
   for (int i=0;i<NDIM*2;i++) {
      cout<<globals.procDomain[i]<<' ';
   }

   // if prompted, init the sorted traverse order
   if (opt) {
      initTraverseOrder();
   }

   // reset all data (neighbours, ghosts, mpi comms etc)
   reset();
}

/*
 * allocate buffers used for mpi communication. 
 * TODO: these buffers take up quite a bit of memory. Some could be combined,
 *       for example the ones used to sync info and the normal ghost buffers
 */
void CdataLL::allocateBuffers() {
   ghostBuffersSend.resize(3);
   pBufferSizes.resize(3);
   gBufferSizes.resize(3);

   /*
    * sending and receiving ghost buffers 
    */
   int totalSize = 0;
   for (Array<CghostData*,NDIM>::iterator i=ghostBuffersSend.begin();i!=ghostBuffersSend.end();i++) {
      vectInt coords = i.position()-1;
      if (globals.procNeighbrs(i.position()) >= 0) {
         *i = new CghostData[calcBufferSize(coords)];
         totalSize += calcBufferSize(coords)*sizeof(CghostData);
         gBufferSizes(i.position()) = calcBufferSize(coords);
      }
   }

   ghostBuffersRecv.resize(3);
   for (Array<CghostData*,NDIM>::iterator i=ghostBuffersRecv.begin();i!=ghostBuffersRecv.end();i++) {
      vectInt coords = i.position()-1;
      if (globals.procNeighbrs(i.position()) >= 0) {
         *i = new CghostData[calcBufferSize(coords)];
         totalSize += calcBufferSize(coords)*sizeof(CghostData);
      }
   }

   /*
    * ghostedParticles contains links to particles that have been ghosted
    * (for later syncing purposes)
    */
   ghostedParticles.resize(3);
   for (Array<Cparticle**,NDIM>::iterator i=ghostedParticles.begin();i!=ghostedParticles.end();i++) {
      vectInt coords = i.position()-1;
      if (globals.procNeighbrs(i.position()) >= 0) {
         *i = new Cparticle*[calcBufferSize(coords)];
         totalSize += calcBufferSize(coords)*sizeof(Cparticle*);
      }
   }

   /*
    * sending and receiving particle buffers 
    */
   particleBuffersSend.resize(3);
   for (Array<Cparticle*,NDIM>::iterator i=particleBuffersSend.begin();i!=particleBuffersSend.end();i++) {
      vectInt coords = i.position()-1;
      if (globals.procNeighbrs(i.position()) >= 0) {
         *i = new Cparticle[calcPBufferSize(coords)];
         totalSize += calcPBufferSize(coords)*sizeof(Cparticle);
         pBufferSizes(i.position()) = calcPBufferSize(coords);
      }
   }

   particleBuffersRecv.resize(3);
   for (Array<Cparticle*,NDIM>::iterator i=particleBuffersRecv.begin();i!=particleBuffersRecv.end();i++) {
      vectInt coords = i.position()-1;
      if (globals.procNeighbrs(i.position()) >= 0) {
         *i = new Cparticle[calcPBufferSize(coords)];
         totalSize += calcPBufferSize(coords)*sizeof(Cparticle);
      }
   }
   
   /*
    * sending and receiving buffers for syncing
    */
   syncBuffersSend.resize(3);
   for (Array<void*,NDIM>::iterator i=syncBuffersSend.begin();i!=syncBuffersSend.end();i++) {
      vectInt coords = i.position()-1;
      if (globals.procNeighbrs(i.position()) >= 0) {
         *i = new unsigned char[sizeof(CghostData)*calcBufferSize(coords)];
         totalSize += calcBufferSize(coords)*sizeof(CghostData);
      }
   }
   syncBuffersRecv.resize(3);
   for (Array<void*,NDIM>::iterator i=syncBuffersRecv.begin();i!=syncBuffersRecv.end();i++) {
      vectInt coords = i.position()-1;
      if (globals.procNeighbrs(i.position()) >= 0) {
         *i = new unsigned char[sizeof(CghostData)*calcBufferSize(coords)];
         totalSize += calcBufferSize(coords)*sizeof(CghostData);
      }
   }
   cout << "size of buffers = "<<totalSize/1024/1024<<" MB"<<endl;
   sendSizesGhosts.resize(3);
   recvSizesGhosts.resize(3);
}


/*
 * reset particle positions, find neighbours, mpi communication to neighbor
 * cpus etc.
 * 
 * user will normally call this whenever particle positions are updated 
 */
void CdataLL::reset() {
   
   tagSortedParticlePointers.clear();

   timeval t1,t2;
   gettimeofday(&t1,NULL);
   int psCapacity = ps.capacity();

   // update domain and sync with neighbouring processes
   updateDomain();

   // check that the sizes of ps, pInfos and pInfoLinks are identical
   if ((ps.size() != pInfos.size())||(ps.size() != pInfoLinks.size())) {
      if (globals.sphStep != 0) {
         cout <<"Error: ps has "<<ps.size()<<" particles, pInfos has "<<pInfos.size()<<" particles, pInfoLinks has "<<pInfoLinks.size()<<" partcles"<<endl;
         exit(-1);
      }
   }

   // ps should not reallocate its memory!
   if (ps.size() > psCapacity) {
      cerr << "Error: dataLL:reset(): ps has grown beyond its capacity! Bad things can happen......I'm outta here!"<<endl;
      exit(-1);
   }
   n = ps.size();

   // recalculate limits
   calcDomainLimits();
   
   // clear bucket list and reinsert particles and ghost particles
   cells.resize(static_cast<vectInt>((rmax-rmin)/(KERNAL_RADIUS*hmax))+8);
#ifdef LIQ_DEM
   dem_cells.resize(static_cast<vectInt>((rmax-rmin)/(2*DEM_RADIUS))+2*3*H/DEM_RADIUS);
   for (Array<vector<Cparticle *>,NDIM>::iterator i=dem_cells.begin();i!=dem_cells.end();i++) {
      i->clear();
   }
#endif

   for (Array<vector<Cparticle *>,NDIM>::iterator i=cells.begin();i!=cells.end();i++) {
      i->clear();
   }
   for (particleContainer::iterator p = ps.begin();p != ps.end();p++) {
      insert(*p);
   }
   for (vector<Cparticle>::iterator p=procGhostParticles.begin();p!= procGhostParticles.end();p++) {
      insert(*p);
   }

   gettimeofday(&t2,NULL);
   globals.wtTotalReset.tv_sec += t2.tv_sec-t1.tv_sec;
   globals.wtTotalReset.tv_usec += t2.tv_usec-t1.tv_usec;
}

void CdataLL::insertNewParticle(Cparticle &newP) {
   ps.push_back(newP);
   n++;
   if (n>MAX_NUM_PARTICLES_PER_CPU) {
      cout<<"insertNewParticle(): number of particles greater than  max limit!!!"<<endl;
      exit(-1);
   }

   /*
    * This is the tricky bit. Make sure that the inserted particle info is also
    *placed in pInfos and pInfoLinks
    */
   Cparticle *p = &(ps.back());
   CpInfo newpInfo;
   newpInfo.p = p;
   newpInfo.key = calcKey(newpInfo.p->r);
   pInfos.push_back(newpInfo);
   CpInfoLink newpInfoLink;
   newpInfoLink.ppInfo = &(pInfos.back());
   pInfoLinks.push_back(newpInfoLink);
   pInfos.back().ppInfoLink = &(pInfoLinks.back());
   insert(*p);
}


void CdataLL::markForDeletion(Cparticle &p) {
   // using a custom tag for deletion. Could be better....
   p.tag = -111;
}

/*
 * deletes all particles with a tag = -111. This could be done better...
 * (delete flag maybe)
 */
void CdataLL::deleteParticles() {
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      Cparticle *thisP = ppInfo->p;
      if (thisP->tag == -111) {
         deleteParticle(ppInfoLink,ppInfo);
         continue;
      }
      ppInfoLink++;
   }
}

/*
 * deletes particle from particle list (ps). Makes sure that the particle
 * infomation is also deleted from pInfoLinks and pInfos
 */
void CdataLL::deleteParticle(vector<CpInfoLink>::iterator &ppInfoLink,CpInfo *ppInfo) {
   Cparticle *thisP = ppInfo->p;
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


/*
 * BIG TODO: This function handles ALL the mpi processing for each timestep.
 * Need to break this up into sub-functions and tidy it all up
 */
void CdataLL::updateDomain() {
   timeval t11,t22;
   gettimeofday(&t11,NULL);

   Array<int,NDIM> ghostIndicies(3);
   Array<int,NDIM> particleIndicies(3);
   ghostIndicies = 0;
   particleIndicies = 0;

   procGhostParticles.clear();

   /*
    * update relevant domain limits (for ghost particles etc) using
    * globals.procDomain.
    */
   vect dmin,dspace,dgmin,dgspace;
#ifdef LIQ_DEM
   vect liq_dgspace,liq_dgmin,dem_dgspace,dem_dgmin;
#endif
   for (int i=0;i<NDIM;i++) {
      dmin[i] = globals.procDomain[i*2];
      dspace[i] = globals.procDomain[i*2+1] - dmin[i];
      vectInt coords = 1;
      coords[i] = 0;
      dgmin[i] = dmin[i] + procGhostMax2H(coords);
      dgspace[i] = dspace[i] - procGhostMax2H(coords);
#ifdef LIQ_DEM
      liq_dgmin[i] = dmin[i] + max(procGhostMax2H(coords),LIQ_DEM_COUPLING_RADIUS);
      dem_dgmin[i] = dmin[i] + 2.0*DEM_RADIUS;
      liq_dgspace[i] = dspace[i] - max(procGhostMax2H(coords),LIQ_DEM_COUPLING_RADIUS);
      dem_dgspace[i] = dspace[i] - 4.0*DEM_RADIUS;
#endif
      coords[i] = 2;
      dgspace[i] -= procGhostMax2H(coords);
#ifdef LIQ_DEM
      liq_dgspace[i] -= max(procGhostMax2H(coords),LIQ_DEM_COUPLING_RADIUS);
#endif
   }
   //cout <<"liq_dgmin = "<<liq_dgmin<<" liq_dgspace = "<<liq_dgspace<<" dem_dgmin = "<<dem_dgmin<<" dem_dgspace = "<<dem_dgspace<<endl;

   /*
    * loop thru particles and add all particle to be sent to
    * particleBuffersSend. Delete these particles from this CPUs list
    */
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;

      ppInfo->neighbrs.clear();
      Cparticle *thisP = ppInfo->p;

      if (thisP->tag == -111) {
         deleteParticle(ppInfoLink,ppInfo);
         continue;
      }
      vectInt outCoords = static_cast<vectInt>((thisP->r-dmin)/dspace+1);
      if (any(outCoords>2)||any(outCoords<0)) {
         cout << "Particle out of range! Tag = "<<thisP->tag<<" iam = "<<thisP->iam<<" position = "<<thisP->r<<" velocityhat = "<<thisP->vhat<<" velocity = "<<thisP->v<<" density = "<<thisP->dens<<" total force = "<<thisP->f<<" pressure force = "<<thisP->fp<<" boundary force = "<<thisP->fb<<endl;
         exit(-1);
      }
      if (globals.procNeighbrs(outCoords) >= 0) {
         //send particle to neighbour
         //add to real particle to send
         if (particleIndicies(outCoords)>=pBufferSizes(outCoords)) {
            cerr << "particleBuffers are full, exiting...."<<endl;
            exit(-1);
         }
         particleBuffersSend(outCoords)[particleIndicies(outCoords)] = *thisP;
         particleIndicies(outCoords)++;
         deleteParticle(ppInfoLink,ppInfo);
         continue;
      }
      ppInfoLink++;
   }

   /*
    * Send these particles to the correct neighbouring cpus
    */
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
     
   /*
    * All particles have arrived. Now loop thru them all and add them to
    * my list. Also, process any special boundary conditions (periodic, ghost etc)
    */
   particleIndicies = particleSizesRecv/sizeof(Cparticle);

   for (Array<int,NDIM>::iterator ap = globals.procNeighbrs.begin();ap != globals.procNeighbrs.end();ap++) {
      vectInt coords = ap.position();
      int neighbr = globals.procNeighbrs(coords);
      if (neighbr >= 0) {
         double maxh = hmax;

         for (int z=0;z<particleIndicies(coords);z++) {
            Cparticle *recvP = &(particleBuffersRecv(coords)[z]);
            ps.push_back(*recvP);
            Cparticle *p = &(ps.back());
            for (int i=0;i<NDIM;i++) {
               if ((PERIODIC[i])&&(coords[i]==0)&&(RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)) {
                  p->r[i] = RMIN[i]-RMAX[i]+p->r[i]; 
                  if (p->iam==sph)
                     p->dens += DENS_DROP[i];
                  #ifdef SLK
                  p->currR[i] = RMIN[i]-RMAX[i]+p->currR[i];
                  #endif
               }
               if ((PERIODIC[i])&&(coords[i]==2)&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP)) {
                  p->r[i] = RMAX[i]-RMIN[i]+p->r[i]; 
                  if (p->iam==sph)
                     p->dens -= DENS_DROP[i];
                  #ifdef SLK
                  p->currR[i] = RMAX[i]-RMIN[i]+p->currR[i];
                  #endif
               }
               if ((GHOST[2*i])&&(coords[i]==0)&&(RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)) {
                  p->r[i] = 2.0*RMIN[i]-p->r[i];
                  if (GHOST[2*i]==1) {
                     p->v[i] = -p->v[i];
                     p->vhat[i] = -p->vhat[i];
                  } else {
                     p->v = -p->v;
                     p->vhat = -p->vhat;
                  }
               }
               if ((GHOST[2*i+1])&&(coords[i]==2)&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP)) {
                  p->r[i] = 2.0*RMAX[i]-p->r[i];
                  if (GHOST[2*i+1]==1) {
                     p->v[i] = -p->v[i];
                     p->vhat[i] = -p->vhat[i];
                  } else {
                     p->v = -p->v;
                     p->vhat = -p->vhat;
                  }
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

   /*
    * now loop thru all the particles again and find all the ghost particles
    * that need to be sent to neighbouring cpus
    */
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();ppInfoLink++) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      ppInfo->neighbrs.clear();
      Cparticle *thisP = ppInfo->p;
#ifdef LIQ_DEM
      vectInt inCoords;
      if (thisP->iam==dem) {
         inCoords = static_cast<vectInt>((thisP->r-dem_dgmin)/dem_dgspace+1);
      } else {
         inCoords = static_cast<vectInt>((thisP->r-liq_dgmin)/liq_dgspace+1);
      }
#else
      vectInt inCoords = static_cast<vectInt>((thisP->r-dgmin)/dgspace+1);
#endif
      for (int i=0;i<NDIM;i++) {
         if (inCoords[i]>2) inCoords[i]=2;
         if (inCoords[i]<0) inCoords[i]=0;
      }
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
               //cout <<"adding type "<<thisP->iam<<" at r = "<<thisP->r<<" to "<<dCoords<<endl; 
               //add to ghost particle array
               if (ghostIndicies(dCoords)>=gBufferSizes(dCoords)) {
                  cerr << "ghostBuffers are full, exiting...."<<endl;
                  cerr << "buffer size is "<<gBufferSizes(dCoords) <<endl;
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

   /*
    * ok. We have found all the ghost and put them in ghostBuffersSend. Now
    * send the data to the correct neighbouring cpus
    */
   Array<int,NDIM> ghostSizesSend(3);
   Array<int,NDIM> ghostSizesRecv(3);
   ghostSizesSend = ghostIndicies*sizeof(CghostData);

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
      
   /*
    * all ghost data is now sent and received. Now loop thru the ghost
    * particles obtained from each neighbour and add them to procGhostParticles, At
    * the same time handle and special boundary conditions.
    */
   ghostIndicies = ghostSizesRecv/sizeof(CghostData);
   recvSizesGhosts = ghostIndicies;

   for (Array<int,NDIM>::iterator ap = globals.procNeighbrs.begin();ap != globals.procNeighbrs.end();ap++) {
      vectInt coords = ap.position();
      int neighbr = globals.procNeighbrs(coords);
      if (neighbr >= 0) {
         //work out ghost max h's
         double maxh = hmax;

         for (int z=0;z<ghostIndicies(coords);z++) {
            Cparticle p;
            p = ghostBuffersRecv(coords)[z];
            for (int i=0;i<NDIM;i++) {
               if ((PERIODIC[i])&&(coords[i]==0)&&(RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)) {
                  p.r[i] = RMIN[i]-RMAX[i]+p.r[i]; 
                  if (p.iam==sph)
                     p.dens += DENS_DROP[i];
                  #ifdef SLK
                  p.currR[i] = RMIN[i]-RMAX[i]+p.currR[i];
                  #endif
               }
               if ((PERIODIC[i])&&(coords[i]==2)&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP)) {
                  p.r[i] = RMAX[i]-RMIN[i]+p.r[i]; 
                  if (p.iam==sph)
                     p.dens -= DENS_DROP[i];
                  #ifdef SLK
                  p.currR[i] = RMAX[i]-RMIN[i]+p.currR[i];
                  #endif
               }
               if ((GHOST[2*i])&&(coords[i]==0)&&(RMIN[i] >= globals.procDomain[i*2]-PSEP)&&(RMIN[i] <= globals.procDomain[i*2]+PSEP)) {
                  p.r[i] = 2.0*RMIN[i]-p.r[i];
                  if (GHOST[2*i]==1) {
                     p.v[i] = -p.v[i];
                     p.vhat[i] = -p.vhat[i];
                  } else {
                     p.v = -p.v;
                     p.vhat = -p.vhat;
                  }
               }
               if ((GHOST[2*i+1])&&(coords[i]==2)&&(RMAX[i] >= globals.procDomain[i*2+1]-PSEP)&&(RMAX[i] <= globals.procDomain[i*2+1]+PSEP)) {
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

         //store max h for this neighbour
         procGhostMax2H(coords) = KERNAL_RADIUS*maxh;
      }
   }

   gettimeofday(&t22,NULL);
   globals.wtTotalUpdateDomain.tv_sec += t22.tv_sec-t11.tv_sec;
   globals.wtTotalUpdateDomain.tv_usec += t22.tv_usec-t11.tv_usec;
}

/*
 * non-blocking mpi send of the particle data in the buffer
 * particleBuffersSend. The function sends the number of particles to be sent before sending the data
 */
void CdataLL::sendRecvParticles(const vectInt coords,const int num,MPI_Request *request,int *particleSizesSend,int *particleSizesRecv) {
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
   MPI_Isend(particleSizesSend,1,MPI_INT,neighbr,300+num,MPI_COMM_WORLD,request);
   MPI_Irecv(particleSizesRecv,1,MPI_INT,neighbr,300+oppNum,MPI_COMM_WORLD,request+1);
   MPI_Isend(particleBuffersSend(coords),*particleSizesSend,MPI_BYTE,neighbr,400+num,MPI_COMM_WORLD,request+2);
   MPI_Irecv(particleBuffersRecv(coords),pBufferSizes(coords)*sizeof(Cparticle),MPI_BYTE,neighbr,400+oppNum,MPI_COMM_WORLD,request+3);
}

/*
 * non-blocking mpi send of the ghost data in the buffer
 * ghostBuffersSend. The function sends the number of particles to be sent before sending the data
 */
void CdataLL::sendRecvGhosts(const vectInt coords,const int num,MPI_Request *request,int *ghostSizesSend,int *ghostSizesRecv) {
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
   MPI_Isend(ghostSizesSend,1,MPI_INT,neighbr,100+num,MPI_COMM_WORLD,request);
   MPI_Irecv(ghostSizesRecv,1,MPI_INT,neighbr,100+oppNum,MPI_COMM_WORLD,request+1);
   MPI_Isend(ghostBuffersSend(coords),*ghostSizesSend,MPI_BYTE,neighbr,200+num,MPI_COMM_WORLD,request+2);
   MPI_Irecv(ghostBuffersRecv(coords),gBufferSizes(coords)*sizeof(CghostData),MPI_BYTE,neighbr,200+oppNum,MPI_COMM_WORLD,request+3);
}

/*
 * for the user-initiated syncing, the number of ghost particles is already
 * know. So no need to send sizes
 */
void CdataLL::sendRecvDataSync(const vectInt coords,const int num,MPI_Request *request,int sendSize, void *sendBuffer, int recvSize, void *recvBuffer) {
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
   MPI_Isend(sendBuffer,sendSize,MPI_BYTE,neighbr,500+num,MPI_COMM_WORLD,request);
   MPI_Irecv(recvBuffer,recvSize,MPI_BYTE,neighbr,500+oppNum,MPI_COMM_WORLD,request+1);
}


void CdataLL::setGlobalTimestep(double dt) {
   double sendbuf = dt;
   double recvbuf[globals.mpiSize];
   MPI_Allgather(&sendbuf, 1, MPI_DOUBLE, recvbuf, 1, MPI_DOUBLE, MPI_COMM_WORLD); 
   globals.dt = dt;
   for (int i=0;i<globals.mpiSize;i++) {
      if (recvbuf[i]<globals.dt) globals.dt = recvbuf[i];
   }
}


/*
 * TODO: should make a general one: "function" over procs
 * TODO: now have two functions for two types (vect and double). 
 * template this to one function
 */
void CdataLL::sumOverProcs(vect *data,int num) {
   
   int newNum = num*sizeof(vect);
   vect *sendbuf = data;
   vect recvbuf[globals.mpiSize*num];

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
}

void CdataLL::sumOverProcs(double *data,int num) {
   
   int newNum = num*sizeof(double);
   double *sendbuf = data;
   double recvbuf[globals.mpiSize*num];

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
}


unsigned long CdataLL::calcKey(vect r) {
   unsigned int rint[NDIM];
   unsigned long key = 0;
   for (int i=0;i<NDIM;i++) {
      rint[i] = (unsigned int)floor((r[i]-rmin[i])*KEY_DIM_MAX/(rmax[i]-rmin[i]));
   }
   unsigned int mask = 0x01<<(KEY_DIM_BITS-1);
   for (int i=0;i<KEY_DIM_BITS;i++) {
      for (int j=0;j<NDIM;j++) {
         key <<= 1;
         key |= (rint[j] & mask) >> (KEY_DIM_BITS-1-i);
      }
      mask >>= 1;
   } 
   return key; 
}

void CdataLL::updateParticles() {
   initPinfos();
   checkPinfos();
}

/*
 * this function clears the particle information contained in pInfos and
 * pInfoLinks and re-generates it from the particle list (ps)
 */
void CdataLL::initPinfos() {
   pInfos.clear();
   pInfoLinks.clear();
   pInfos.reserve(ps.capacity());
   pInfoLinks.reserve(ps.capacity());
   cout << "size of pInfos and pInfoLinks = "<<(pInfos.capacity()*sizeof(CpInfo) + pInfoLinks.capacity()*sizeof(CpInfoLink))/1024/1024<<" MB"<<endl;
   for (particleContainer::iterator thisP = ps.begin();thisP!=ps.end();thisP++) {
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

/*
 * Sort pInfos using insertion sort. Make sure that the sorting
 * maintains consistency with pInfoLink
 */
void CdataLL::sortPInfos() {
   for (vector<CpInfo>::iterator i=pInfos.begin()+1;i!=pInfos.end();i++) {
      for (vector<CpInfo>::iterator j=pInfos.begin();j!=i;j++) {
         if (j->key > i->key) {
            CpInfo tmpInfo = *j;
            *j = *i;
            j->ppInfoLink->ppInfo = &(*j);
            vector<CpInfo>::iterator jplus1 = j+1;
            for (vector<CpInfo>::iterator k=i;k!=jplus1;k--) {
               *k = *(k-1);
               k->ppInfoLink->ppInfo = &(*k);
            }
            *jplus1 = tmpInfo;
            jplus1->ppInfoLink->ppInfo = &(*jplus1);
         }
      }
   }
}

/*
 * for each member of pInfoLinks, the particle information pointed to by ppInfo
 * should have a class member ppInfoLink that points back to the original entry in pInfoLinks
 */
void CdataLL::checkPinfos() {
   particleContainer::iterator p = ps.begin();
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();ppInfoLink++) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      if (ppInfo->ppInfoLink != &(*ppInfoLink)) {
         cerr <<"Error: CdataLL::initTraverseOrder: ppInfo and ppInfoLink are not consistant!"<<endl;
         exit(-1);
      }
      if (ppInfo->p != &(*p)) {
         cerr <<"Error: CdataLL::initTraverseOrder: ppInfo, ppInfoLink and p are not consistant!"<<endl;
         exit(-1);
      }
      p++;
   }
}


void CdataLL::initTraverseOrder() {
   initPinfos();   
   sortPInfos();
   checkPinfos();
}

void CdataLL::calcTraverseOrder() {
   particleContainer::iterator p = ps.begin();
   for (vector<CpInfoLink>::iterator ppInfoLink = pInfoLinks.begin();ppInfoLink!=pInfoLinks.end();ppInfoLink++) {
      CpInfo *ppInfo = ppInfoLink->ppInfo;
      if (ppInfo->ppInfoLink != &(*ppInfoLink)) {
         cerr <<"Error: CdataLL::calcTraverseOrder: ppInfo and ppInfoLink are not consistant!"<<endl;
         exit(-1);
      }
      if (ppInfo->p != &(*p)) {
         cerr <<"Error: CdataLL::calcTraverseOrder: ppInfo, ppInfoLink and p are not consistant!"<<endl;
         exit(-1);
      }
      ppInfo->key = calcKey(ppInfo->p->r);
      ppInfo->neighbrs.clear();
      p++;
   }
   sortPInfos();
}

void CdataLL::calcDomainLimits() {
   for (int i=0;i<NDIM;i++) {
      rmin[i] = globals.procDomain[i*2];
      rmax[i] = globals.procDomain[i*2+1];
   }
   hmax = -100000;
   for (particleContainer::iterator thisP = ps.begin();thisP!=ps.end();thisP++) {
      for (int j=0;j<NDIM;j++) {
         if (thisP->r[j] < rmin[j]) { rmin[j] = thisP->r[j]; }
         if (thisP->r[j] > rmax[j]) { rmax[j] = thisP->r[j]; }
      }
     
      if (thisP->h > hmax) hmax = thisP->h;
   }
   for (int i=0;i<NDIM;i++) {
      vectInt coords = 1;
      coords[i] = 0;
      if (procGhostMax2H(coords)<hmax) procGhostMax2H(coords)=KERNAL_RADIUS*hmax;
      coords[i] = 2;
      if (procGhostMax2H(coords)<hmax) procGhostMax2H(coords)=KERNAL_RADIUS*hmax;
   }
}


void CdataLL::insert(Cparticle &p) {
#ifdef LIQ_DEM
   if (p.iam==dem) {
      vectInt cellI = dem_getCellI(p.r);
      dem_cells(cellI).push_back(&p);
   } else {
      vectInt cellI = getCellI(p.r);
      cells(cellI).push_back(&p);
   }
#else
   vectInt cellI = getCellI(p.r);
   cells(cellI).push_back(&p);
#endif
}

void CdataLL::addIfNeighbr(Cparticle &_p,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd) {
   for (int i=0;i<_toAdd.size();i++) {
      if (_p.is_neighbr(*(_toAdd[i]))) {
         _neighbrs.push_back(_toAdd[i]);
      }
   }
}

void CdataLL::addIfInRadius(Cparticle &_p,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd,const double radius2) {
   for (int i=0;i<_toAdd.size();i++) {
      if (len2(_p.r-_toAdd[i]->r)<radius2) {
         _neighbrs.push_back(_toAdd[i]);
      }
   }
}

void CdataLL::addAll(Cparticle &_p,vector<Cparticle *> &_neighbrs,vector<Cparticle *> &_toAdd) {
   for (int i=0;i<_toAdd.size();i++) {
      _neighbrs.push_back(_toAdd[i]);
   }
}

/*
 * finds bucket of given particle. Loops over neighbouring buckets in order to
 * find neighbouring particles.
 */
void CdataLL::calcNeighbours(vector<Cparticle *> &_neighbrs,Cparticle &_p) {
#ifdef LIQ_DEM
   if (_p.iam==dem) {
   const int num_cells = 1;
   vectInt cellI = dem_getCellI(_p.r);
   vectInt tcellI = cellI;
   for (int i=cellI[0]-num_cells;i<=cellI[0]+num_cells;i++) {
      tcellI[0] = i;
      if (NDIM > 1) {
         for (int j=cellI[1]-num_cells;j<=cellI[1]+num_cells;j++) {
            tcellI[1] = j;
            if (NDIM > 2) {
               for (int k=cellI[2]-num_cells;k<=cellI[2]+num_cells;k++) {
                  tcellI[2] = k;
                  addIfNeighbr(_p,_neighbrs,dem_cells(tcellI));
               }
            } else {
               addIfNeighbr(_p,_neighbrs,dem_cells(tcellI));
            }
         }
      } else {
         addIfNeighbr(_p,_neighbrs,dem_cells(tcellI));
      } 
   }
   } else {
#endif
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
#ifdef LIQ_DEM
   }
#endif
}


void CdataLL::calcNeighboursAtRadius(vector<Cparticle *> &_neighbrs,Cparticle &_p,const double radius,const double radius2) {
   const int num_cells = int(ceil(0.5*radius/hmax));
   vectInt cellI = getCellI(_p.r);
   vectInt tcellI = cellI;
   for (int i=cellI[0]-num_cells;i<=cellI[0]+num_cells;i++) {
      tcellI[0] = i;
      if (NDIM > 1) {
         for (int j=cellI[1]-num_cells;j<=cellI[1]+num_cells;j++) {
            tcellI[1] = j;
            if (NDIM > 2) {
               for (int k=cellI[2]-num_cells;k<=cellI[2]+num_cells;k++) {
                  tcellI[2] = k;
                  addIfInRadius(_p,_neighbrs,cells(tcellI),radius2);
               }
            } else {
               addIfInRadius(_p,_neighbrs,cells(tcellI),radius2);
            }
         }
      } else {
         addIfInRadius(_p,_neighbrs,cells(tcellI),radius);
      } 
   }
}


#ifdef LIQ_DEM
void CdataLL::calcNeighboursCoupling(vector<Cparticle *> &_neighbrs,Cparticle &_p) {
   const int num_cells = int(ceil(0.5*LIQ_DEM_COUPLING_RADIUS/hmax));
   vectInt cellI = getCellI(_p.r);
   vectInt tcellI = cellI;
   for (int i=cellI[0]-num_cells;i<=cellI[0]+num_cells;i++) {
      tcellI[0] = i;
      if (NDIM > 1) {
         for (int j=cellI[1]-num_cells;j<=cellI[1]+num_cells;j++) {
            tcellI[1] = j;
            if (NDIM > 2) {
               for (int k=cellI[2]-num_cells;k<=cellI[2]+num_cells;k++) {
                  tcellI[2] = k;
                  addIfInRadius(_p,_neighbrs,cells(tcellI),pow(LIQ_DEM_COUPLING_RADIUS,2));
               }
            } else {
               addIfInRadius(_p,_neighbrs,cells(tcellI),pow(LIQ_DEM_COUPLING_RADIUS,2));
            }
         }
      } else {
         addIfInRadius(_p,_neighbrs,cells(tcellI),pow(LIQ_DEM_COUPLING_RADIUS,2));
      } 
   }
   if (_p.iam!=dem) {
      const int num_cells = int(ceil(0.5*LIQ_DEM_COUPLING_RADIUS/DEM_RADIUS));
      vectInt cellI = dem_getCellI(_p.r);
      vectInt tcellI = cellI;
      for (int i=cellI[0]-num_cells;i<=cellI[0]+num_cells;i++) {
         tcellI[0] = i;
         const int iPlus2 = pow(i-cellI[0]+1,2);
         const int iMinus2 = pow(i-cellI[0]-1,2);
         if (NDIM > 1) {
            for (int j=cellI[1]-num_cells;j<=cellI[1]+num_cells;j++) {
               const int jPlus2 = pow(j+1-cellI[1],2);
               const int jMinus2 = pow(j-1-cellI[1],2);
               tcellI[1] = j;
               if (NDIM > 2) {
                  for (int k=cellI[2]-num_cells;k<=cellI[2]+num_cells;k++) {
                     tcellI[2] = k;
                     const int rPlus2 = iPlus2+jPlus2+pow(k-cellI[2]+1,2);
                     const int rMinus2 = iMinus2+jMinus2+pow(k-cellI[2]-1,2);
                     if (rPlus2<pow(LIQ_DEM_COUPLING_RADIUS/DEM_RADIUS,2)) {
                        addAll(_p,_neighbrs,dem_cells(tcellI));
                     } else if (rMinus2<pow(LIQ_DEM_COUPLING_RADIUS/DEM_RADIUS,2)) {
                        addIfInRadius(_p,_neighbrs,dem_cells(tcellI),pow(LIQ_DEM_COUPLING_RADIUS,2));
                     }
                  }
               } else {
                  const int rPlus2 = iPlus2+jPlus2;
                  const int rMinus2 = iMinus2+jMinus2;
                  if (rPlus2<pow(LIQ_DEM_COUPLING_RADIUS/DEM_RADIUS,2)) {
                     addAll(_p,_neighbrs,dem_cells(tcellI));
                  } else if (rMinus2<pow(LIQ_DEM_COUPLING_RADIUS/DEM_RADIUS,2)) {
                     addIfInRadius(_p,_neighbrs,dem_cells(tcellI),pow(LIQ_DEM_COUPLING_RADIUS,2));
                  }
               }
            }
         } else {
            addIfInRadius(_p,_neighbrs,dem_cells(tcellI),pow(LIQ_DEM_COUPLING_RADIUS,2));
         } 
      }
   }
}
#endif

Cparticle *CdataLL::getParticleFromTag(int tag) {
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
      cerr << "Error in CdataLL::getParticleFromTag(): the tag aint here!!!"<<endl;
      exit(-1);
   }
   return tagSortedParticlePointers[tag-1].p;
}

CpInfo *CdataLL::getPinfoFromTag(int tag) {
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
      cerr << "Error in CdataLL::getPinfoFromTag(): the tag aint here!!!"<<endl;
      exit(-1);
   }
   return tagSortedPinfoPointers[tag-1].pInfo;
}

int CdataLL::calcBufferSize(vectInt coords) {
   int bufferSize = 1.0;
   int bufferSizeDem = 1.0;
#ifdef HACK_TO_SAVE_MEM
   if (coords[2]!=-1) {
      const double nnxx = pow((1-POROSITY)/DEM_VOL,1.0/NDIM);
      const double demSep = 1.0/nnxx;
      bufferSizeDem *= pow(2.0*DEM_RADIUS,NDIM)/pow(demSep,NDIM);
   }
#endif
   for (int i=0;i<NDIM;i++) {
      if (coords[i]==0) {
#ifdef LIQ_DEM
         bufferSizeDem *= ceil((RMAX[i]-RMIN[i])/(2.0*DEM_RADIUS));
#endif
         bufferSize *= ceil((RMAX[i]-RMIN[i])/PSEP);
      } else {
#ifdef LIQ_DEM
         bufferSizeDem *= 2;
         bufferSize *= ceil(LIQ_DEM_COUPLING_RADIUS/PSEP);
#endif
         bufferSize *= ceil(2.0*H/PSEP);
      }
   }
#ifdef LIQ_DEM
   return 1.5*(bufferSize+min(bufferSizeDem,MAX_NUM_DEM_PARTICLES));
#else
   return 1.5*bufferSize;
#endif
}    
int CdataLL::calcPBufferSize(vectInt coords) {
   int bufferSize = 1;
   int bufferSizeDem = 1;
#ifdef HACK_TO_SAVE_MEM
   if (coords[2]!=-1) {
      const double nnxx = pow((1-POROSITY)/DEM_VOL,1.0/NDIM);
      const double demSep = 1.0/nnxx;
      bufferSizeDem *= pow(2.0*DEM_RADIUS,NDIM)/pow(demSep,NDIM);
   }
#endif
   for (int i=0;i<NDIM;i++) {
      if (coords[i]==0) {
#ifdef LIQ_DEM
         bufferSizeDem *= ceil((RMAX[i]-RMIN[i])/(2.0*DEM_RADIUS));
#endif
         bufferSize *= ceil((RMAX[i]-RMIN[i])/PSEP);
      } else {
#ifdef LIQ_DEM
         bufferSizeDem *= 3;
#endif
         bufferSize *= 3;
      }
   }
#ifdef LIQ_DEM
   return 1.5*(bufferSize+min(bufferSizeDem,MAX_NUM_DEM_PARTICLES));
#else
   return 1.5*bufferSize;
#endif
}   
         
   
   
