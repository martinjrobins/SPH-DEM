
#include "misc.h"
#include "misc.impl.h"
#include "dataLL.impl.h"

#include "sph.h"
#include "io_data_vtk.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



const double PAIR_DBL_R_MULT2 = pow(1.2,2);
const double PAIR_DBL_MAX_TIME =50;
const int PAIR_DBL_SAMPLES_PER_OUTSTEP = 10000;

using namespace Nsph;


#ifdef LIQ_DEM

void Nmisc::addRandomDEMparticlesCyl(CdataLL *data,const vect origin,const double radius,const double height,int n) {
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd1);
   
   //gsl_rng_set(rng,1);
   gsl_rng_set(rng,time(NULL));
   
   cout << "will put in "<<n<<" dem particles..."<<endl;

   vect minRange = origin;
   vect maxRange = origin;
   minRange[1] -= radius;
   minRange[0] -= radius;
   maxRange[1] += radius;
   maxRange[0] += radius;
   maxRange[2] += height;

   Cparticle p;
   const double radius2 = pow(radius,2);
   int origNumParticles = data->getParticles()->size();
   while(data->getParticles()->size() < n+origNumParticles) {
      for (int ii=0;ii<NDIM;ii++) {
         p.r[ii] = gsl_ran_flat(rng,minRange[ii],maxRange[ii]);
      }
      if (pow(p.r[0],2)+pow(p.r[1],2) > radius2) continue;
      p.tag = data->getParticles()->size()+1;
      p.dens = DEM_DENS;
      p.h = DEM_RADIUS;
      p.mass = DEM_VOL*DEM_DENS;
      p.v = 0.0;
      p.vhat = p.v;
      p.iam = dem;
      vector<Cparticle *> neighbrs;
      data->calcNeighbours(neighbrs,p);
      bool good = true;
      for (vector<Cparticle *>::iterator ii=neighbrs.begin();ii!=neighbrs.end();ii++) {
         if (((*ii)->iam == dem)&&(len2((*ii)->r-p.r)<=4.0*pow(DEM_RADIUS,2))) {
            good = false;
         }
      }
      if (good) {
         cout <<".";
         data->insertNewParticle(p);
         data->reset();
      }
   }
   cout<<endl;
   gsl_rng_free(rng);
}

void Nmisc::addRandomDEMparticles(CdataLL *data,const double minRange[NDIM],const double maxRange[NDIM],const double porosity) {
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd1);
   
   //gsl_rng_set(rng,1);
   gsl_rng_set(rng,time(NULL));
   
   double volume = 1;
   for (int i=0;i<NDIM;i++) {
      volume *= maxRange[i]-minRange[i];
   }
   
   const int n = (int)((1-porosity)*volume/((4.0/3.0)*PI*pow(DEM_RADIUS,3)));
   cout << "will put in "<<n<<" dem particles..."<<endl;

   Cparticle p;
   int origNumParticles = data->getParticles()->size();
   while(data->getParticles()->size() < n+origNumParticles) {
      for (int ii=0;ii<NDIM;ii++) {
         p.r[ii] = gsl_ran_flat(rng,minRange[ii],maxRange[ii]);
      }
      p.tag = data->getParticles()->size()+1;
      p.dens = DEM_DENS;
      p.h = DEM_RADIUS;
      p.mass = DEM_VOL*DEM_DENS;
      p.v = 0.0;
      p.vhat = p.v;
      p.iam = dem;
      vector<Cparticle *> neighbrs;
      data->calcNeighbours(neighbrs,p);
      bool good = true;
      for (vector<Cparticle *>::iterator ii=neighbrs.begin();ii!=neighbrs.end();ii++) {
         if (((*ii)->iam == dem)&&(len2((*ii)->r-p.r)<=4.0*pow(DEM_RADIUS,2))) {
            good = false;
         }
      }
      if (good) {
         cout <<".";
         data->insertNewParticle(p);
         data->reset();
      }
   }
   cout<<endl;
   gsl_rng_free(rng);
}

void Nmisc::addGridDEMparticles(vector<Cparticle> &ps,const double minRange[NDIM],const double maxRange[NDIM],const double porosity) {
  
   double volume = 1;
   for (int i=0;i<NDIM;i++) {
      volume *= maxRange[i]-minRange[i];
   }
   
   const double n = ((1-porosity)*volume/DEM_VOL);
   cout << "will put in "<<n<<" dem particles..."<<endl;
   const double dperUnit = pow(n/volume,1.0/NDIM);
   double dsep = 1.0/dperUnit;
   if (dsep<2.0*DEM_RADIUS) {
      cout << "Porosity is too low to fit in the particles!!!"<<endl;
      cout <<"dsep = "<<dsep<<" dperUnit = "<<dperUnit<<endl;
      exit(-1);
   }
   const int dnx = dperUnit*(maxRange[0]-minRange[0]);
   const int dny = dperUnit*(maxRange[1]-minRange[1]);
#ifdef _3D_
   const int dnz = dperUnit*(maxRange[2]-minRange[2]);
#endif
   dsep = (maxRange[1]-minRange[1])/dny;

   Cparticle p;
   for (int i=0;i<dnx;i++) {
      for (int j=0;j<dny;j++) {
#ifdef _3D_
         for (int k=0;k<dnz;k++) {
            p.r = (i+0.5)*dsep+minRange[0],(j+0.5)*dsep+minRange[1],(k+0.5)*dsep+minRange[2];
#else

            p.r = (i+0.5)*dsep+minRange[0],(j+0.5)*dsep+minRange[1];
#endif
            p.tag = ps.size()+1;
            p.dens = DEM_DENS;
            p.h = DEM_RADIUS;
            p.mass = DEM_VOL*DEM_DENS;
            p.v = 0.0;
            p.vhat = p.v;
            p.iam = dem;
            cout <<".";
            ps.push_back(p);
#ifdef _3D_
         }
#endif
      }
   }
   cout<<endl;
}

#endif



void Nmisc::boundaryLine(vector<Cparticle> &ps,vect *points,vect *normals,bool *concave,int numPoints,enum iamTypes iam) {
   Cparticle p;
   for (int i = 0;i < numPoints-1;i++) {
      double lineLength = len(points[i+1]-points[i]);
      vect m = (points[i+1]-points[i])/lineLength;
      vect c = points[i];
      int numPointsInLine = int(floor(lineLength/PSEP+0.5));
      double thisPSEP = lineLength/numPointsInLine;
      if (i>0) {
         p.concave = concave[i];
         p.norm2 = normals[i-1];
      }
      for (int j=0;j<numPointsInLine;j++) {
         p.r = j*thisPSEP*m + c;
         p.tag = ps.size()+1;
         p.dens = DENS;
         p.mass = pow(BFAC*PSEP,NDIM)*DENS;
         p.h = H;
         p.v = 0,0;
         p.iam = iam;
         p.norm1 = normals[i];
         p.dist = 0;
         ps.push_back(p);
         p.norm2 = 0.0;
      }
   }
}

#ifdef _3D_
void Nmisc::boundaryPlane(vector<Cparticle> &ps,vect side1,vect side2,vect neighbrNorm1,vect neighbrNorm2,vect neighbrNorm3,vect neighbrNorm4,vect initPos,double thisPsep,bool flipNorm) {
   double len1 = len(side1);
   double len2 = len(side2);
   int wN1 = floor(len1/thisPsep);
   double wPSEP1 = len1/wN1;
   int wN2 = floor(len2/thisPsep);
   double wPSEP2 = len2/wN2;
   vect norm1 = side1/len1;
   vect norm2 = side2/len2;
   Cparticle p;

   int startI,startJ,endI,endJ;

   if (all(neighbrNorm1==0)) {
      startI = 1;
   } else {
      startI = 0;
   }
   if (all(neighbrNorm2==0)) {
      startJ = 1;
   } else {
      startJ = 0;
   }
   if (all(neighbrNorm3==0)) {
      endI = wN1-1;
   } else {
      endI = wN1;
   }
   if (all(neighbrNorm4==0)) {
      endJ = wN2-1;
   } else {
      endJ = wN2;
   }
      
   for (int i=startI;i<=endI;i++) {
      for (int j=startJ;j<=endJ;j++) {
         p.tag = ps.size()+1;
         p.r = initPos;
         p.r += i*wPSEP1*norm1;
         p.r += j*wPSEP2*norm2;
         p.dens = DENS;
         p.mass = pow(PSEP,NDIM)*DENS;
         p.h = HFAC*pow(p.mass/p.dens,1.0/NDIM);
         p.v = 0.0;
         p.vhat = p.v;
         p.iam = boundary;
         int count = 1;
         p.norm1 = cross(norm1,norm2);
         if (flipNorm) p.norm1 = -p.norm1;
         p.norm3 = norm2;
         if (i==0) {
            vect normtmp = (p.norm1*count + neighbrNorm1)/(count+1);
            p.norm3 = cross(p.norm1,neighbrNorm1);
            p.norm1 = normtmp;
         }
         if (i==wN1) {
            vect normtmp = (p.norm1*count + neighbrNorm3)/(count+1);
            p.norm3 = cross(p.norm1,neighbrNorm3);
            p.norm1 = normtmp;
         }
         if (j==0) {
            vect normtmp = (p.norm1*count + neighbrNorm2)/(count+1);
            p.norm3 = cross(p.norm1,neighbrNorm2);
            p.norm1 = normtmp;
         }
         if (j==wN2) {
            vect normtmp = (p.norm1*count + neighbrNorm4)/(count+1);
            p.norm3 = cross(p.norm1,neighbrNorm4);
            p.norm1 = normtmp;
         }
         p.norm1 = p.norm1/len(p.norm1);
         p.norm3 = p.norm3/len(p.norm3);
         p.norm2 = 0;
         ps.push_back(p);
      }
   }
}

void Nmisc::boundaryCircle(vector<Cparticle> &ps,const vect &origin,const double radiusMin,const double radiusMax,const vect normal) {
   Cparticle p;
   const double dRadius = radiusMax-radiusMin;
   const double NR = ceil(dRadius/(BFAC*PSEP));
   const double RPSEP = dRadius/NR;
   const vect n(1,0,0);
   for (int i=0;i<=NR;i++) {
      const double r = i*RPSEP+radiusMin;
      double nc = ceil(2*PI*r/(BFAC*PSEP));
      if (nc == 0) nc++;
      const double TPSEP = 2.0*PI/nc;
      for (int j=0;j<nc;j++) {
         const double theta = j*TPSEP;
         vect newN;
         newN[0] = n[0]*cos(theta)-n[1]*sin(theta);
         newN[1] = n[0]*sin(theta)+n[1]*cos(theta);
         newN[2] = 0;
         p.r = r*newN + origin;
         p.tag = ps.size()+1;
         p.dens = DENS;
         p.mass = pow(BFAC*PSEP,NDIM)*DENS;
         p.h = H;
         p.v = 0,0,0;
         p.iam = boundary;
         p.norm1 = normal;
#ifdef _3D_
         p.norm3 = cross(newN,p.norm1);
#endif
         if (i==NR) {
            p.norm2[0] = -cos(theta);
            p.norm2[1] = -sin(theta);
            p.norm2[2] = 0;
            p.concave = 1;
         } else if ((i==0)&&(radiusMin>0)) {
            p.norm2[0] = -cos(theta);
            p.norm2[1] = -sin(theta);
            p.norm2[2] = 0;
            p.concave = 1;
         } else {
            p.norm2 = 0,0,0;
         }
         p.dist = 0;
         ps.push_back(p);
      }
   }
}

void Nmisc::boundaryCylinderNoTopBottom(vector<Cparticle> &ps,const vect &origin,const double radius,const double height) {
   Cparticle p;
   const double NH = ceil(height/(BFAC*PSEP));
   const double HPSEP = height/NH;
   const vect n(1,0,0);
   const vect nt(0,0,1);
   double nc = ceil(2*PI*radius/(BFAC*PSEP));
   if (nc == 0) nc++;
   const double TPSEP = 2.0*PI/nc;
   for (int i=0;i<NH;i++) {
      for (int j=0;j<nc;j++) {
         const double theta = j*TPSEP;
         vect newN;
         newN[0] = n[0]*cos(theta)-n[1]*sin(theta);
         newN[1] = n[0]*sin(theta)+n[1]*cos(theta);
         p.r = radius*newN + nt*i*HPSEP + origin;
         p.tag = ps.size()+1;
         p.dens = DENS;
         p.mass = pow(BFAC*PSEP,NDIM)*DENS;
         p.h = H;
         p.v = 0,0,0;
         p.iam = boundary;
         p.norm1[0] = -cos(theta);
         p.norm1[1] = -sin(theta);
         p.norm1[2] = 0;
         p.norm3[0] = sin(theta);
         p.norm3[1] = -cos(theta);
         p.norm3[2] = 0;
         p.norm2 = 0,0,0;
         p.dist = 0;
         ps.push_back(p);
      }
   }
}


void Nmisc::boundaryCylinder(vector<Cparticle> &ps,const vect &origin,const double radius,const double height) {
  vect normal = 0.0;
  normal[2] = 1.0;
  boundaryCircle(ps,origin,0,radius,normal);
  vect newOrigin = origin;
  const double newHeight = height - BFAC*PSEP;
  newOrigin[2] += BFAC*PSEP;
  cout << "creating cylinder with origin "<<newOrigin<<endl;
  boundaryCylinderNoTopBottom(ps,newOrigin,radius,newHeight);
}
void Nmisc::sphCylinder(vector<Cparticle> &ps,const vect &origin) {
  Cparticle p;
  const vect nt(0,0,1);
  const vect n(1,0,0);
  for (int k = 0;k<NZ;k++) {
     for (int j = 0;j<NX;j++) {
        const double radius = j*PSEP;
        int ntheta = ceil(2.0*PI*radius/PSEP);
        if (ntheta==0) ntheta++;
        const double dtheta = 2.0*PI/ntheta;
        for (int i=0;i<ntheta;i++) {
         const double theta = i*dtheta;
         vect newN;
         newN[0] = n[0]*cos(theta)-n[1]*sin(theta);
         newN[1] = n[0]*sin(theta)+n[1]*cos(theta);
         newN[2] = 0;
         p.r = radius*newN + nt*k*PSEP + origin;
         p.tag = ps.size()+1;
         p.dens = DENS;
         p.mass = pow(PSEP,NDIM)*DENS;
         p.h = H;
         p.v = 0,0,0;
         p.iam = sph;
         ps.push_back(p);
        }
     }
  }
}
#endif


inline bool ifInteriorSph(Cparticle &p) {
   return (p.iam == sph)&&(abs(p.r[0])<RMAX[0]*0.75)&&(abs(p.r[1])<RMAX[1]*0.75);
}
void Nmisc::pairDoubling(CdataLL *data,Sgi::slist<CpairData> &pairData,vector<ofstream *> &openFiles,gsl_rng *rng,Cio_data_vtk *io) {
   //pair doubling stuff
   //go through existing pairs
   Sgi::slist<CpairData>::iterator i = pairData.begin();
   Sgi::slist<CpairData>::iterator i_before;
   while (i != pairData.end()) {
      Cparticle *p1 = data->getParticleFromTag(i->p1);
      Cparticle *p2 = data->getParticleFromTag(i->p2);
      if ((p1->tag!=i->p1)||(p2->tag!=i->p2)) {
         cerr << "Error in pairDoubling: getParticleFromTag didn't work!"<<endl;
         exit(-1);
      }
      double dt = data->globals.time-i->initT;
      double r2 = len2(p1->r-p2->r);
      if ((r2>i->finalR2)||(dt>PAIR_DBL_MAX_TIME)) {
         double lyap;
         if (dt>PAIR_DBL_MAX_TIME) {
            lyap = 0;
         } else {
            lyap = (0.5/dt) * log(r2/i->initR2);
            *(i->outFile) << sqrt(i->initR2) <<" "<< lyap <<" "<<sqrt(r2)<<" "<<dt<<endl;
         }
         if (i==pairData.begin()) {
            i_before = i;
            i++;
            pairData.erase(i_before);
         } else {
            i++;
            pairData.erase_after(i_before);
         }
      } else {
         i_before = i;
         i++;
      }
   }

   //open new output file
   string outputFN = "PairDoublingTimes";
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",data->globals.outstep);
   outputFN = io->getFilename()+outputFN+strTimestep+".dat";
   ofstream *newFile = new ofstream(outputFN.c_str());
   *newFile <<"#initial_dist lyaponov_exponent final_dist time_taken"<<endl;
   openFiles.push_back(newFile);
   
   //init new pairs
   vector<Cparticle> *ips = data->getParticles();
   Cparticle *p1;
   Cparticle *p2;
   for (int i=0;i<PAIR_DBL_SAMPLES_PER_OUTSTEP;i++) {
      getRandomPair<ifInteriorSph>(*ips,rng,&p1,&p2);
      
      CpairData newPD;
      newPD.p1 = p1->tag;
      newPD.p2 = p2->tag;
      newPD.outFile = newFile;
      newPD.initT = data->globals.time;
      newPD.initR2 = len2(p1->r-p2->r);
      newPD.finalR2 = newPD.initR2*PAIR_DBL_R_MULT2;
      pairData.push_front(newPD);
   }

   if (data->globals.mpiRank==0) cout <<"\t\tThere are now "<<pairData.size()<<" pairs being tracked."<<endl;
}
 
inline void calcParticleLyap(Cparticle &p, vector<Cparticle *> &neighbrs,CglobalVars &g,CdataLL &d2) {
   Array<double,2> coeff(NDIM,NDIM);

   double chisq;

   int n = neighbrs.size();

   gsl_matrix *X = gsl_matrix_alloc(n,NDIM);
   gsl_vector *y = gsl_vector_alloc(n);

   gsl_vector *c = gsl_vector_alloc(NDIM);
   gsl_matrix *cov = gsl_matrix_alloc(NDIM,NDIM);
            
   gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n,NDIM);

   vector<vect> dx2(n);
   for (int i=0;i<n;i++) {
      Cparticle *pn = neighbrs[i];
      vect dx = pn->r-p.r;
      Cparticle *pn2 = d2.getParticleFromTag(pn->tag);
      for (int j=0;j<NDIM;j++) {
         if (PERIODIC[j]) {
            if (abs(pn2->r[j]-p.r[j]+RMAX[j]-RMIN[j])<abs(pn2->r[j]-p.r[j])) {
               (dx2[i])[j] = pn2->r[j]-p.r[j]+RMAX[j]-RMIN[j];
            } else if (abs(pn2->r[j]-p.r[j]-RMAX[j]+RMIN[j])<abs(pn2->r[j]-p.r[j])) {
               (dx2[i])[j] = pn2->r[j]-p.r[j]-RMAX[j]+RMIN[j];
            } else {
               (dx2[i])[j] = pn2->r[j]-p.r[j];
            }
         } else {
            (dx2[i])[j] = pn2->r[j]-p.r[j];
         }
         gsl_matrix_set(X,i,j,dx[j]);
      }
   }

   for (int j=0;j<NDIM;j++) {
      for (int i=0;i<n;i++) {
         gsl_vector_set(y,i,(dx2[i])[j]);
      }
            
      gsl_multifit_linear(X,y,c,cov,&chisq,work);

      for (int i=0;i<NDIM;i++) {
         coeff(i,j) = gsl_vector_get(c,i);
      }
   }
   gsl_multifit_linear_free(work);
   gsl_vector_free(c);
   gsl_vector_free(y);
   gsl_matrix_free(X);
   gsl_matrix_free(cov);

   //TODO: this is 2d!, change to 3d!
   // multiply coeff by its transverse
   double aa,bb,cc,dd;
   aa = pow(coeff(0,0),2) + pow(coeff(1,0),2);
   bb = coeff(0,0)*coeff(0,1) + coeff(1,0)*coeff(1,1);
   cc = bb;
   dd = pow(coeff(0,1),2) + pow(coeff(1,1),2);

   //calculate eigenvalues of resultant matrix
   double eigen[2];
   eigen[0] = (aa+dd + sqrt(pow(aa+dd,2) - 4.0*(aa*dd-bb*cc)))/2.0;
   eigen[1] = (aa+dd - sqrt(pow(aa+dd,2) - 4.0*(aa*dd-bb*cc)))/2.0;
   
   double DLE = 1.0/(2.0*abs(d2.globals.time - g.time))*log(max(eigen[0],eigen[1]));
   
   p.tmp = DLE;
}
           


void Nmisc::calcLyap(CdataLL *d1,CdataLL *d2,Cio_data_vtk *io,double dt) {
   vector<double> ftle,bftle;
   ftle.resize(d1->getParticles()->size());
   bftle.resize(d2->getParticles()->size());

   d1->neighboursGroup<CdataLL,calcParticleLyap,ifSph>(*d2);
   particleContainer *p1 = d1->getParticles();
   int n = p1->size();
   for (int i=0;i<n;i++) {
      ftle[i] = (*p1)[i].tmp;
   }

   /*
   d2->neighboursGroup<CdataLL,calcParticleLyap,ifSph>(*d1);
   particleContainer *p2 = d2->getParticles();
   n = p2->size();
   for (int i=0;i<n;i++) {
      bftle[i] = (*p2)[i].tmp;
   }
   */

   stringstream ss;
   ss << int(abs(dt));
   string outName;
   if (dt<0) {
      outName = "bftle_DT"+ss.str()+"_";
   } else {
      outName = "ftle_DT"+ss.str()+"_";
   }
   io->writeAux(d1->globals.outstep,*(d1->getParticles()),ftle,outName.c_str(),&(d1->globals));
   //io->writeAux(d2->globals.outstep,*(d2->getParticles()),bftle,"bftle",&(d2->globals));
}


void Nmisc::calcPoincare(CdataLL *d1,CdataLL *d2,Cio_data_vtk *io,int numParticles,int numPeriods) {
   particleContainer ps;
   Cparticle p;
   int tags[numParticles];
   Cparticle *tmpPs[numParticles];

   gsl_rng *rng;
   rng = gsl_rng_alloc(gsl_rng_ranlxd1);
   //gsl_rng_set(rng,1);
   gsl_rng_set(rng,time(NULL));

   vector<Cparticle> *p1 = d1->getParticles();
   int totalNumParticles = p1->size();

   for (int i=0;i<numParticles;i++) {
      int np = (int)floor(gsl_rng_uniform(rng)*totalNumParticles);
      while (!(Nsph::ifSph((*p1)[np]))) {
         np = (int)floor(gsl_rng_uniform(rng)*totalNumParticles);
      }
      tmpPs[i] = &((*p1)[np]);
   }

   for (int i=0;i<numPeriods;i++) {
     for (int j=0;j<numParticles;j++) { 
       Cparticle *p2 = d2->getParticleFromTag(tmpPs[j]->tag);
       ps.push_back(*p2);

       vector<Cparticle *> neighbrs;
       d1->calcNeighbours(neighbrs,*p2);
       double closestR2 = pow((RMAX[0]-RMIN[0]),2);
       Cparticle *closestP = neighbrs[0];
       for (vector<Cparticle *>::iterator pNeighbr = neighbrs.begin();pNeighbr!=neighbrs.end();pNeighbr++) {
          if (Nsph::ifSph(**pNeighbr)) {
            double r2 = len2((*pNeighbr)->r-p2->r);
            if (r2<closestR2) {
               closestP = (*pNeighbr);
               closestR2 = r2;
            }
          }
       }
       tmpPs[j] = closestP; 
     }
   }
   io->writeAuxNoData(d1->globals.outstep,ps,"poincare",&(d1->globals));
}
   
   

int Nmisc::coordsToNum(vectInt coords,vectInt split) {
   int num = 0;
   for (int j=NDIM-1;j>=0;j--) {
      num = num*split[j] + coords[j];
   }
   return num;
}

vectInt Nmisc::numToCoords(int num,vectInt split) {
   int tmp = num;
   vect coords;
   for (int j=0;j<NDIM;j++) {
      coords[j] = tmp%split[j];
      tmp /= split[j];
   }
   return coords; 
}

void Nmisc::setupEntireDomain(double procDomain[NDIM*2],Array<int,NDIM> &procNeighbrs) {
   vector<vector<double> > vprocDomain(1);
   vector<Array<int,NDIM> > vprocNeighbrs(1);
   vector<particleContainer > vps;
   vectInt split;
   for (int i=0;i<NDIM;i++) {
      split[i] = 1;
   }
   particleContainer pps;
   Nmisc::splitDomain(pps,split,vps,vprocDomain,vprocNeighbrs);
   for (int i=0;i<NDIM*2;i++) {
      procDomain[i] = vprocDomain[0][i];
   }
   procNeighbrs = vprocNeighbrs[0];
}

void Nmisc::divideParticlesToCPUs(particleContainer &ps,vector<particleContainer > &vps,vector<vector<double> > &vprocDomain) {
   const int ncpu = vprocDomain.size();
   const int n = ps.size();
   for (int i=0;i<ncpu;i++) {
      cout << "for cpu "<<i<<" vprocdomain = ";
      for (int b=0;b<3*2;b++) cout<<vprocDomain[i][b]<<' ';
      cout<<endl;
      for (int j=0;j<n;j++) {
         bool in = true;
         for (int k=0;k<NDIM;k++) {
            if ((ps[j].r[k] <= vprocDomain[i][2*k])||(ps[j].r[k] > vprocDomain[i][2*k+1])) {
               in = false;
            }
         }
         if (in) {
            //cout <<"putting "<<ps[j].r<<" into cpu "<<i<<endl;
            vps[i].push_back(ps[j]);
         }
      }
   }
}

void Nmisc::splitDomain(particleContainer &ps,vectInt &split,vector<particleContainer > &vps,vector<vector<double> > &vprocDomain,vector<Array<int,NDIM> > &vprocNeighbrs) {
 
   int nProc = product(split);
   //cout << "we have "<<nProc<<" processors"<<endl;
   vps.resize(nProc);
   vprocDomain.resize(nProc);
   vprocNeighbrs.resize(nProc);


   vect spacing;

   for (int i=0;i<NDIM;i++) {
      spacing[i] = (RMAX[i]-RMIN[i])/split[i];
   }
   //cout <<"spacings are: "<<spacing[0]<<' '<<spacing[1]<<endl;

   Array<int,NDIM> neighbrMap(split+2);
   for (Array<int,NDIM>::iterator p = neighbrMap.begin();p!=neighbrMap.end();p++) {
      vectInt coords = p.position()-1;
      for (int i=0;i<NDIM;i++) {
         if (coords[i] < 0) {
            if (PERIODIC[i]) {
               coords[i] = split[i]-1;
            } else if (GHOST[2*i]) {
               coords[i] = 0;
            } else {
               *p = -1;
               break;
            }
         } else if (coords[i] >= split[i]) {
            if (PERIODIC[i]) {
               coords[i] = 0;
            } else if (GHOST[2*i+1]) {
               coords[i] = split[i]-1;
            } else {
               *p = -1;
               //cout << "coords = "<<coords<<endl;
               break;
            }
         }
         *p = coordsToNum(coords,split);
      }
   }
   //cout << "Neighbour map is "<<neighbrMap<<endl;

   for (int i=0;i<nProc;i++) {
      vprocDomain[i].resize(NDIM*2);
      vprocNeighbrs[i].resize(3);

      vectInt coords = numToCoords(i,split);
      //cout <<"doing coords = "<<coords<<endl;
      for (int j=0;j<NDIM;j++) {
         if ((coords[j]==0)&&(!PERIODIC[j])&&(!GHOST[2*j])) {
            //vprocDomain[i][2*j] = RMIN[j]-2*KERNAL_RADIUS*H;
            vprocDomain[i][2*j] = RMIN[j]-0.25*(RMAX[j]-RMIN[j]);
         } else {
            vprocDomain[i][2*j] = spacing[j]*coords[j]+RMIN[j];
         }
         if ((coords[j]==split[j]-1)&&(!PERIODIC[j])&&(!GHOST[2*j+1])) {
            //vprocDomain[i][2*j+1] = RMAX[j]+2*KERNAL_RADIUS*H;
            vprocDomain[i][2*j+1] = RMAX[j]+0.25*(RMAX[j]-RMIN[j]);
         } else {
            vprocDomain[i][2*j+1] = spacing[j]*(coords[j]+1)+RMIN[j];
         }
      }
      RectDomain<NDIM> subdomain(coords,coords+2);
      vprocNeighbrs[i] = neighbrMap(subdomain);
#ifdef _2D_
      vprocNeighbrs[i](1,1) = -1;
#else
      vprocNeighbrs[i](1,1,1) = -1;
#endif
      //cout <<"neighbrs = "<<vprocNeighbrs[i]<<endl;
      //cout <<"so n(0,1) = "<<vprocNeighbrs[i](0,1)<<endl;
   }

   for (particleContainer::iterator thisP=ps.begin();thisP!=ps.end();thisP++) {
      vectInt coords;
      for (int j=0;j<NDIM;j++) {
         coords[j] = int((thisP->r[j]-RMIN[j])/spacing[j]);
         if (coords[j]>= split[j]) coords[j] = split[j]-1;
         if (coords[j]< 0) coords[j] = 0;
      }
      vps[coordsToNum(coords,split)].push_back(*thisP);
   }
}
