
#ifndef MISC_IMPL_H
#define MISC_IMPL_H

#include "misc.h"

const int VEL_STRUCT_SAMPLES_PER_OUTSTEP = 1000;
const int POS_STRUCT_SAMPLES_PER_OUTSTEP = 1000;
const int FULL_VEL_STRUCT_RES = 3000;

template<bool ifFunct(Cparticle &p)>
void Nmisc::getRandomNeighbrPair(CdataLL *data,gsl_rng *rng,Cparticle **p1,Cparticle **p2){
   
   vector<Cparticle> *ips = data->getParticles();
   int numParticles = ips->size();
   int ip1 = (int)floor(gsl_rng_uniform(rng)*numParticles);
   while (!ifFunct((*ips)[ip1])) {
      ip1 = (int)floor(gsl_rng_uniform(rng)*numParticles);
   }
   vector<Cparticle *> neighbrs;
   data->calcNeighbours(neighbrs,(*ips)[ip1]);
   int numNeighbrs = neighbrs.size();
   int ip2 = (int)floor(gsl_rng_uniform(rng)*numNeighbrs);
   int j=0;
   while (!ifFunct(*(neighbrs[ip2]))) {
      ip2 = (int)floor(gsl_rng_uniform(rng)*numNeighbrs);
      j++;
   }
   *p1 = &((*ips)[ip1]);
   *p2 = neighbrs[ip2]; 
}



template<bool ifFunct(Cparticle &p)>
void Nmisc::getRandomPair(particleContainer &ps,gsl_rng *rng,Cparticle **p1,Cparticle **p2){
   int finished = 0;
   while (!finished) {
      int numParticles = ps.size();
      int ip1 = (int)floor(gsl_rng_uniform(rng)*numParticles);
      while (!ifFunct(ps[ip1])) {
         ip1 = (int)floor(gsl_rng_uniform(rng)*numParticles);
      }
      const int numBins = int(floor(log((RMAX[0]-RMIN[0])/PSEP)/log(2.0)));
      const int bin = (int)floor(gsl_rng_uniform(rng)*numBins);
      const double r2 = pow((RMAX[0]-RMIN[0])/pow(2.0,numBins-bin),2);
      const double r2_floor = pow((RMAX[0]-RMIN[0])/pow(2.0,numBins-bin+0.5),2);
      const double r2_ceil = pow((RMAX[0]-RMIN[0])/pow(2.0,numBins-bin-0.5),2);
      vector<Cparticle *> candidates;
      for (int j=0;j<numParticles;j++) {
         double thisR2 = len2(ps[j].r-ps[ip1].r);
         if ((thisR2>=r2_floor)&&(thisR2<=r2_ceil)&&ifFunct(ps[j])) {
            candidates.push_back(&(ps[j]));
         }
      }
      int ip2;   
      if (!candidates.empty()) {
         ip2 = (int)floor(gsl_rng_uniform(rng)*candidates.size());
         finished = 1;
      } else {
         cout <<"No candidates!: r = "<<ps[ip1].r<<" r_floor = "<<sqrt(r2_floor)<<" r2_ceil = "<<sqrt(r2_ceil)<<endl;
         continue;
      }
      *p1 = &(ps[ip1]);
      *p2 = candidates[ip2]; 
   }
}

inline void calcVelocityMode(Cparticle &p,CglobalVars &g,vectInt &coords) {
   vect refr=2.0*p.r/(RMAX[1]-RMIN[1]);
   double lenCentre = len(refr);
   refr = (refr+1.0)/2.0;
   //double v2m2rd = 1.0/(RMAX[0]-RMIN[0])*len2(p.v)*pow(p.mass,2)/p.dens;
   double v2m2rd;
   if (PERIODIC[0]) {
      v2m2rd = len2(p.v)*pow(p.mass,2)/((RMAX[1]-RMIN[1])*p.dens);
   } else {
      if (lenCentre<1) {
         v2m2rd = 0.5*(1.0-cos(PI*(lenCentre+1)))*(0.5*p.mass*len2(p.v)-(g.eKE/g.nSph))*p.mass/(0.5*(RMAX[1]-RMIN[1])*p.dens);
      } else {
         v2m2rd = 0.0;
      }
   }
   //v2m2rd = len2(p.v)*pow(p.mass,2)/((RMAX[1]-RMIN[1])*p.dens);
   //double v2m2rd = pow(p.mass*p.v[1],2)/((RMAX[1]-RMIN[1])*p.dens);
   double loc = 2.0*PI*dot(coords,refr);
   //cout <<"spectrum with v = "<<p.v<<"mass = "<<p.mass<<" dens = "<<p.dens<<" refr = "<<refr<<endl;
   g.tmpArray[0] += v2m2rd*cos(loc); 
   g.tmpArray[1] += v2m2rd*sin(loc);
}

template<bool ifFunct(Cparticle &p)>
void Nmisc::calcVelocityStructure(CdataLL *data,gsl_rng *rng,Cio_data_vtk *io) {
   string filename = "VelStruct";
   //open new output file
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   string outputFN = io->getFilename()+filename+strTimestep+".dat";
   ofstream newFile(outputFN.c_str(),ios::out|ios::trunc);
   if (data->globals.mpiRank==0) cout <<"calculating velocity structure info for neighbours. Putting into file: "<<outputFN<<endl;
   newFile <<"#dvx dvy drx dry"<<endl;
   
   Cparticle *p1;
   Cparticle *p2;
   for (int i=0;i<VEL_STRUCT_SAMPLES_PER_OUTSTEP;i++) {
      getRandomNeighbrPair<ifFunct>(data,rng,&p1,&p2);
      //getRandomPair<ifFunct>(*ips,rng,&p1,&p2);
      vect dr = p1->r-p2->r;
      vect dv = p1->v-p2->v;
      newFile <<dv[0]<<' '<<dv[1]<<' '<<dr[0]<<' '<<dr[1]<<endl;
   }
   newFile.close();
}

template<bool ifFunct(Cparticle &p)>
void Nmisc::calcFullVelocityStructure(CdataLL *data,Cio_data_vtk *io) {
   string filename = "VelStructFull";
   //open new output file
   if (data->globals.mpiRank==0) cout <<"calculating full velocity structure info"<<endl;
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   string outputFN = io->getFilename()+filename+strTimestep+".dat";
   ofstream newFile(outputFN.c_str(),ios::out|ios::trunc);
   newFile <<"#r v v2 v3 v12"<<endl;
   
   vector<double> resv(FULL_VEL_STRUCT_RES,0);
   vector<double> resv2(FULL_VEL_STRUCT_RES,0);
   vector<double> resv3(FULL_VEL_STRUCT_RES,0);
   vector<double> resv12(FULL_VEL_STRUCT_RES,0);
   vector<unsigned int> counter(FULL_VEL_STRUCT_RES,0);


   vector<Cparticle> *ips = data->getParticles();
   int numParticles = ips->size();

   for (int i=0;i<numParticles;i++) {
      Cparticle *p1 = &((*ips)[i]);
      if (!ifFunct(*p1)) continue;
      for (int j=0;j<numParticles;j++) {
         Cparticle *p2 = &((*ips)[j]);
         if (!ifFunct(*p2)) continue;
         
         vect drVect = p1->r-p2->r; 
         double dr = len(drVect);
         double dv = abs(dot(p1->v-p2->v,drVect)/dr);

         if (dr<(RMAX[1]-RMIN[1])) {
            double dv2 = pow(dv,2);
            double dv3 = dv2*dv;
            double dv12 = dv3*pow(dv2,2);
            int index = int(dr*(FULL_VEL_STRUCT_RES/(RMAX[1]-RMIN[1]))+0.5);
            counter[index]++;
            resv[index] += dv;
            resv2[index] += dv2;
            resv3[index] += dv3;
            resv12[index] += dv12;
         }
      }
   }

   if (data->globals.mpiRank==0) cout <<"writing full velocity structure info to file: "<<outputFN<<endl;

   for (int i=0;i<FULL_VEL_STRUCT_RES;i++) {
      if (counter[i]>0) {
         newFile <<i*((RMAX[1]-RMIN[1])/FULL_VEL_STRUCT_RES)<<' '<<resv[i]/counter[i]<<' '<<resv2[i]/counter[i]<<' '<<resv3[i]/counter[i]<<' '<<resv12[i]/counter[i]<<endl;
      }
   }
   newFile.close();
}

    
template<bool ifFunct(Cparticle &p)>
void Nmisc::calcPositionStructure(CdataLL *data,gsl_rng *rng,Cio_data_vtk *io) {
   //open new output file
   string filename = "PosStruct";
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   string outputFN = io->getFilename()+filename+strTimestep+".dat";
   ofstream newFile(outputFN.c_str(),ios::out|ios::trunc);
   newFile <<"# drx dry"<<endl;

   Cparticle *p1;
   Cparticle *p2;
   for (int i=0;i<POS_STRUCT_SAMPLES_PER_OUTSTEP;i++) {
      getRandomNeighbrPair<ifFunct>(data,rng,&p1,&p2);
      
      vect dr = p2->r-p1->r;
      //TODO: this is 2D specific!! 
      newFile <<dr[0]<<' '<<dr[1]<<endl;
   }
   newFile.close();
}



template<bool ifFunct(Cparticle &p)>
Array<myComplex,NDIM> *Nmisc::calcVelocitySpectrum(CdataLL *data, const int numModes) {
   const int numModes2 = int(pow(double(numModes),2));
   Array<myComplex,NDIM> *modes = new Array<myComplex,NDIM>(numModes);
   *modes = 0;
   for (Array<myComplex,NDIM>::iterator i=modes->begin();i!=modes->end();i++) {
      vectInt coords = i.position();
      int mag2 = int(len2(coords));
      if (mag2<numModes2) {
         data->globals.tmpArray[0] = 0;
         data->globals.tmpArray[1] = 0;
         data->traverse<vectInt,calcVelocityMode,ifFunct>(coords);
         *i = myComplex(data->globals.tmpArray[0],data->globals.tmpArray[1]);
         //*i = sqrt(pow(data->globals.tmpArray[0],2) + pow(data->globals.tmpArray[1],2));
         //*i = sqrt(pow(data->globals.tmpArray[0],2));
      }
   }
   return modes;
}


template<bool ifFunct(Cparticle &p)>
void Nmisc::outputVelocitySpectrum(CdataLL *data,Cio_data_vtk *io) {
   const int numModes = int(0.4*(RMAX[1]-RMIN[1])/PSEP);

   Array<myComplex,NDIM> *modes = calcVelocitySpectrum<ifFunct>(data,numModes); 
   
   vector<double> oneD(numModes,0);
   vector<int> nums(numModes,0);
   for (Array<myComplex,NDIM>::iterator i=modes->begin();i!=modes->end();i++) {
      vectInt coords = i.position();
      int mag = int(len(coords));
      if (mag<numModes) {
         nums[mag] += 1;
         oneD[mag] += abs(*i);
      }
   }
   vector<int>::iterator inum = nums.begin();
   for (vector<double>::iterator i=oneD.begin();i!=oneD.end();i++) {
      *i /= *inum; 
      inum++;
   }
   
   ofstream fo;
   string outputFN = "VelSpec";
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   string outputFN_1D = io->getFilename()+outputFN+"1D"+strTimestep+".dat";
   string outputFN_2D = io->getFilename()+outputFN+"2D"+strTimestep+".dat";

   cout <<"writing 1D spectrum to file: "<<outputFN_1D<<endl;
   fo.open(outputFN_1D.c_str(),ios::out|ios::trunc);
   fo << "# mode 1DSpectrum"<<endl;
   int n = oneD.size();
   for (int i=0;i<n;i++) {
      fo<<i<<' '<<oneD[i]<<endl;
   }
   fo.close();

#ifdef _2D_
   cout <<"writing 2D spectrum to file: "<<outputFN_2D<<endl;
   fo.open(outputFN_2D.c_str(),ios::out|ios::trunc);
   fo << "# "<<numModes<<" x "<<numModes<<" matrix"<<endl;
   for (int i=0;i<numModes;i++) {
      for (int j=0;j<numModes;j++) {
         fo<<' '<<real((*modes)(i,j));
      }
      fo<<endl;
   }
   for (int i=0;i<numModes;i++) {
      for (int j=0;j<numModes;j++) {
         fo<<' '<<imag((*modes)(i,j));
      }
      fo<<endl;
   }
   fo.close();
#endif

   delete modes;
}

template<double specFunct(Cparticle &p,CglobalVars &g)>
inline void calcMode(Cparticle &p,CglobalVars &g,vectInt &coords) {
   vect refr=2.0*p.r/(RMAX[1]-RMIN[1]);
   double lenCentre = len(refr);
   refr = (refr+1.0)/2.0;
   //double v2m2rd = 1.0/(RMAX[0]-RMIN[0])*len2(p.v)*pow(p.mass,2)/p.dens;
   double v2m2rd;
   if (PERIODIC[0]) {
      v2m2rd = 2.0*specFunct(p,g)*p.mass/((RMAX[1]-RMIN[1])*p.dens);
      //v2m2rd = len2(p.v)*pow(p.mass,2)/((RMAX[1]-RMIN[1])*p.dens);
   } else {
      if (lenCentre<1) {
         v2m2rd = (1.0-cos(PI*(lenCentre+1)))*specFunct(p,g)*p.mass/((RMAX[1]-RMIN[1])*p.dens);
         //v2m2rd = 0.5*(1.0-cos(PI*(lenCentre+1)))*(0.5*p.mass*len2(p.v)-(g.eKE/g.nSph))*p.mass/(0.5*(RMAX[1]-RMIN[1])*p.dens);
      } else {
         v2m2rd = 0.0;
      }
   }
   //v2m2rd = len2(p.v)*pow(p.mass,2)/((RMAX[1]-RMIN[1])*p.dens);
   //double v2m2rd = pow(p.mass*p.v[1],2)/((RMAX[1]-RMIN[1])*p.dens);
   double loc = 2.0*PI*dot(coords,refr);
   //cout <<p.tag<<" ";
   //cout <<"spectrum with v = "<<p.v<<"mass = "<<p.mass<<" dens = "<<p.dens<<" refr = "<<refr<<endl;
   g.tmpArray[0] += v2m2rd*cos(loc); 
   g.tmpArray[1] += v2m2rd*sin(loc);
}


template<bool ifFunct(Cparticle &p),double specFunct(Cparticle &p,CglobalVars &g)>
Array<myComplex,NDIM> *Nmisc::calcSpectrum(CdataLL *data, const int numModes) {
   vectInt numModesVect = numModes;
   for (int i=1;i<NDIM;i++) {
      numModesVect[i] = numModesVect[i]*2-1;
   }
   vectInt lbounds = 0;
   for (int i=1;i<NDIM;i++) {
      lbounds[i] = -numModes+1;
   }
   const int numModes2 = int(pow(double(numModes),2));
   Array<myComplex,NDIM> *modes = new Array<myComplex,NDIM>(lbounds,numModesVect);
   *modes = 0;
   for (Array<myComplex,NDIM>::iterator i=modes->begin();i!=modes->end();i++) {
      vectInt coords = i.position();
      int mag2 = int(len2(coords));
      if (mag2<numModes2) {
         //cout << "calculating mode = "<<coords<<endl;
         data->globals.tmpArray[0] = 0;
         data->globals.tmpArray[1] = 0;
         data->traverse<vectInt,calcMode<specFunct>,ifFunct>(coords);
         *i = myComplex(data->globals.tmpArray[0],data->globals.tmpArray[1]);
         //*i = sqrt(pow(data->globals.tmpArray[0],2) + pow(data->globals.tmpArray[1],2));
         //*i = sqrt(pow(data->globals.tmpArray[0],2));
      }
   }
   return modes;
}



template<bool ifFunct(Cparticle &p),void specFunct(Cparticle &p,CglobalVars &g,vectInt &coords)>
void Nmisc::setSpectrum(CdataLL *data, const int numModes) {
   //const vectInt numModesVect = numModes*2-1;
   //const vectInt lbounds = -numModes+1;
   const vectInt numModesVect = numModes;
   const vectInt lbounds = 0;
   const int numModes2 = int(pow(double(numModes),2));
   Array<myComplex,NDIM> *modes = new Array<myComplex,NDIM>(lbounds,numModesVect);

   for (Array<myComplex,NDIM>::iterator i=modes->begin();i!=modes->end();i++) {
      vectInt coords = i.position();
      int mag2 = int(len2(coords));
      if (mag2<numModes2) {
         //cout << "setting mode = "<<coords<<endl;
         data->traverse<vectInt,specFunct,ifFunct>(coords);
      }
   }
}


template<bool ifFunct(Cparticle &p),double specFunct(Cparticle &p,CglobalVars &g)>
void Nmisc::outputSpectrum(CdataLL *data,Cio_data_vtk *io,string name) {
   const int numModes = int(0.4*(RMAX[1]-RMIN[1])/PSEP);

   Array<myComplex,NDIM> *modes = calcSpectrum<ifFunct,specFunct>(data,numModes); 
   
   vector<double> oneD(numModes,0);
   vector<int> nums(numModes,0);
   for (Array<myComplex,NDIM>::iterator i=modes->begin();i!=modes->end();i++) {
      vectInt coords = i.position();
      int mag = int(len(coords));
      if (mag<numModes) {
         nums[mag] += 1;
         oneD[mag] += abs(*i);
      }
   }
   vector<int>::iterator inum = nums.begin();
   for (vector<double>::iterator i=oneD.begin();i!=oneD.end();i++) {
      *i /= *inum; 
      inum++;
   }
   
   ofstream fo;
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   string outputFN_1D = io->getFilename()+name+"1D"+strTimestep+".dat";
   string outputFN_2D = io->getFilename()+name+"2D"+strTimestep+".dat";

   cout <<"writing 1D spectrum to file: "<<outputFN_1D<<endl;
   fo.open(outputFN_1D.c_str(),ios::out|ios::trunc);
   fo << "# mode 1DSpectrum"<<endl;
   int n = oneD.size();
   for (int i=0;i<n;i++) {
      fo<<i<<' '<<oneD[i]<<endl;
   }
   fo.close();

#ifdef _2D_
   cout <<"writing 2D spectrum to file: "<<outputFN_2D<<endl;
   fo.open(outputFN_2D.c_str(),ios::out|ios::trunc);
   fo << "# "<<numModes<<" x "<<numModes<<" matrix"<<endl;
   for (int i=0;i<numModes;i++) {
      for (int j=-numModes+1;j<numModes;j++) {
         fo<<' '<<real((*modes)(i,j));
      }
      fo<<endl;
   }
   for (int i=0;i<numModes;i++) {
      for (int j=-numModes+1;j<numModes;j++) {
         fo<<' '<<imag((*modes)(i,j));
      }
      fo<<endl;
   }
   fo.close();
#endif

   delete modes;
}
#endif
