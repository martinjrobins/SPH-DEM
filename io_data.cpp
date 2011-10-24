#include "io_data.h"
    
vector<Cparticle> *Cio_data::read(int timestep) {
   H5PartFile *reader = H5PartOpenFile(filename.c_str(),H5PART_READ);
   if (H5PartGetNumSteps(reader) < timestep) {
      cerr << "There is not "<<timestep<<"timesteps in the file"<<endl;
      H5PartCloseFile(reader);
      exit(-1);
   }
   H5PartSetStep(reader,timestep);
   int n = H5PartGetNumParticles(reader);
   
   vector<Cparticle> *output = new vector<Cparticle>(n);

   double r[NDIM][n],v[NDIM][n],mass[n],u[n];
   
   for (int i=0;i<NDIM;i++) {
      if (!H5PartReadDataFloat64(reader,dimNames[i],r[i])) {
         cerr << "Error reading particle position data" << endl;
         exit(-1);
      }
   }
   for (int i=0;i<NDIM;i++) {
      if (!H5PartReadDataFloat64(reader,vdimNames[i],v[i])) {
         cerr << "Error reading particle velocity data" << endl;
         exit(-1);
      }
   }
   if (!(H5PartReadDataFloat64(reader,"mass",mass))&&
        (H5PartReadDataFloat64(reader,"u",u))) {
      cerr << "Error reading particle data in file" << endl;
      exit(-1);
   }
   for (int i=0;i<n;i++) {
      for (int j=0;j<NDIM;i++) {
         (*output)[i].r[j] = r[j][i];
         (*output)[i].v[j] = r[j][i];
      }
      (*output)[i].mass = mass[i];
      (*output)[i].u = u[i];
   }
   H5PartCloseFile(reader);
   return output;
}

void Cio_data::write(vector<Cparticle> &ps, int timestep) {
   H5PartFile *writer = H5PartOpenFile(filename.c_str(),H5PART_WRITE);

   H5PartSetStep(writer,timestep);
   int n = ps.size();
   H5PartSetNumParticles(writer,n); 

   double r[NDIM][n],v[NDIM][n],mass[n],u[n];
   long long tag[n];
   
   double zero[n];
   for (int i=0;i<n;i++) {
      zero[i] = 0.0;
      for (int j=0;j<NDIM;j++) {
         r[j][i] = ps[i].r[j];
         v[j][i] = ps[i].v[j];
      }
      mass[i] = ps[i].mass;
      u[i] = ps[i].u;
      tag[i] = ps[i].tag;
   }

   for (int i=0;i<NDIM;i++) {
      H5PartWriteDataFloat64(writer,dimNames[i],r[i]);
      H5PartWriteDataFloat64(writer,vdimNames[i],v[i]);
   }
   for (int i=NDIM;i<3;i++) {
      H5PartWriteDataFloat64(writer,dimNames[i],zero);
      H5PartWriteDataFloat64(writer,vdimNames[i],zero);
   }
   H5PartWriteDataFloat64(writer,"mass",mass);
   H5PartWriteDataFloat64(writer,"u",u);
   H5PartWriteDataInt64(writer,"id",tag);
}
