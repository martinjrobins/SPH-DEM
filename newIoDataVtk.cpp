#include "io_data_vtk.h"

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkLongArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <cstdlib>
#include <boost/tuple/tuple.hpp>


template <class double>
class typeMap {
   typedef vtkDoubleArray vtkArrayType
   typedef double vtkTupleType 
   static int components = 1;
}

template <class vect>
class typeMap {
   typedef vtkDoubleArray vtkArrayType
   typedef double[NDIM] vtkTupleType
   static int components = NDIM;
}

template <class int>
class typeMap {
   typedef vtkIntArray vtkArrayType
   typedef int vtkTupleType 
   static int components = 1;
}

template <class particleT>
inline void readProcessTuple(int i, vtkUnstructuredGrid *vtkData, vector<particleT>::iterator p, const null_type&) {};

template <class particleT, class H, class T>
inline void readProcessTuple(int i, vtkUnstructuredGrid *vtkData, vector<particleT>::iterator p, cons<H,T> vars) {
   typedef typeMap<vars.get_head()::type>::vtkArrayType vtkType
   vtkType *vtkArray = static_cast<vtkType *> (dataset->GetPointData()->GetArray(vars.get_head().name));
   vars.get_head()(*p) = vtkArray->getTuple(i);
   readRestartProcessTuple(i,vtkData,data,vars.get_tail());
}

template <class particleT>
inline void writeProcessTuple(int i, vtkUnstructuredGrid *vtkData, vector<particleT>::iterator p, const null_type&) {};

template <class particleT, class H, class T>
inline void writeProcessTuple(int i, vtkUnstructuredGrid *vtkData, vector<particleT>::iterator p, cons<H,T> vars) {
   typedef typeMap<H::type>::vtkArrayType vtkType
   typedef typeMap<H::type>::vtkTupleType tupleType
   vtkType *vtkArray = static_cast<vtkType *> (vtkData->GetPointData()->GetArray(vars.get_head().name));
   typleType tuple = vars.get_head()(*p);
   vtkArray->SetTuple(i, &tuple);
   writeProcessTuple(i,vtkData,data,vars.get_tail());
}

template <class particleT>
inline void writeInitProcessTuple(int n, vtkUnstructuredGrid *vtkData, const null_type&) {};

template <class particleT, class H, class T>
inline void writeInitProcessTuple(int n, vtkUnstructuredGrid *vtkData, cons<H,T> vars) {
   typedef typeMap<H::type>::vtkArrayType vtkType
   typedef typeMap<H::type>::vtkTupleType tupleType
   const int components = typeMap<H::type>::components;
   vtkType *vtkArray = vtkType::New();
   vtkArray->SetNumberOfComponents(components);
   vtkArray->SetNumberOfTuples(n);
   vtkArray->SetName(H::name);
   vtkData->GetPointData()->AddArray(vtkArray);
   writeInitProcessTuple(i,vtkData,vars.get_tail());
}

vtkXMLUnstructuredGridReader *openReader(int timestep,string filename,CglobalVars *globals) {
   vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   char strRank[20];
   sprintf(strRank,"%d",globals->mpiRank);
   string theFilename = filename+strTimestep+"_"+strRank+".vtu";
   cout << "Reading file " << theFilename << endl;
   reader->SetFileName(theFilename.c_str());
   reader->Update(); 
   return reader;
}

vtkXMLUnstructuredGridWriter *openWriter(int timestep,string filename,CglobalVars *globals) {
   vtkXMLPUnstructuredGridWriter *writer = vtkXMLPUnstructuredGridWriter::New();
   writer->SetNumberOfPieces(globals->mpiSize);
   writer->SetStartPiece(globals->mpiRank);
   writer->SetEndPiece(globals->mpiRank);
   writer->SetDataModeToBinary();
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   string outFilename = filename+strTimestep+".pvtu";
   writer->SetFileName(outFilename.c_str());
   //writer->SetTimeStep(timestep);
   return writer;
}


template <class particleT,class varsT>
void readParticles(string filename, int timestep,vector<particleT> *outPs,CglobalVars *globals,varsT vars) {
   vtkXMLUnstructuredGridReader *reader = openFile(timestep,filename);
   vtkUnstructuredGrid *vtkData = reader->GetOutput();
   
   int n = dataset->GetNumberOfPoints();
   if (globals->mpiRank==0) cout << "Found "<<n<<" points in the file..." << endl;

   outPs->resize(n);
   
   vector<particleT>::iterator psIter = outPs->begin();
   for (int i=0;i<n;i++) {
      readProcessTuple(i, *vtkData, psIter, vars);
      psIter++;
   }
   reader->Delete();
}

void deleteArrays(vtkPoints *pts) {
  const n = pts->GetNumberOfArrays();
  for (int i=0;i<n;i++) {
     pts->GetArray(i)->Delete();
  }
}

template <class particleT,class varsT>
void writeParticles(string filename, int timestep,vector<particleT> *outPs,CglobalVars *globals,varsT vars) {
   vtkUnstructuredGrid *vtkData);
   
   const int n = outPs->size();

   vtkPoints *newPts = vtkPoints::New();
   vtkCellArray *newCells = vtkCellArray::New();
   newPts->SetNumberOfPoints(n);
   vtkUnstructuredGrid *vtkData = vtkUnstructuredGrid::New();
   vtkData->SetPoints(newPts);
   vtkData->SetCells(1,newCells);

   writeInitProcessTuple(n, vtkUnstructuredGrid *vtkData, cons<H,T> vars) {
   
   vector<particleT>::iterator psIter = outPs->begin();
   for (int i=0;i<n;i++) {
      newCells->InsertNextCell(1);
      newCells->InsertCellPoint(i);
      double r = psIter->r;
      newPts->SetPoint(i,&r);
      writeProcessTuple(i, vtkData, psIter, vars);
      psIter++;
   }

   if (globals->mpiRank==0) cout << "\tWriting "<<n<<" particles to file: "<<outFilename<<endl;

   vtkXMLUnstructuredGridWriter writer = openWriter(timestep,filename,globals);
   writer->SetInput(vtkData);
   writer->Write();

   deleteArrays(newPts);

   writer->Delete();
   dataset->Delete();
   newPts->Delete();
   newCells->Delete();
}

template<class particleT>
void Cio_data_vtk::readRestartDem(int timestep,vector<particleT> *outPs,CglobalVars *globals) {
   outPs->reserve(MAX_NUM_PARTICLES_PER_CPU);
   readParticles(filename+"Restart",timestep,outPs,globals,particleT::restartVars());
}

template<class particleT>
void Cio_data_vtk::writeRestart(int timestep,vector<particleT> *outPs,CglobalVars *globals) {
   writeParticles(filename+"Restart",timestep,outPs,globals,particleT::restartVars());
}

void Cio_data_vtk::readRestartSph(int timestep,vector<sphParticle> *outPs,CglobalVars *globals) {
    
   
   vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   char strRank[20];
   sprintf(strRank,"%d",globals->mpiRank);
   string theFilename = filename+"Restart"+strTimestep+"_"+strRank+".vtu";
   cout << "Reading file " << theFilename << endl;
   reader->SetFileName(theFilename.c_str());
   reader->Update(); 
   vtkUnstructuredGrid *vtkData = reader->GetOutput();
   //reader->SetOutput(dataset);
   
   int n = dataset->GetNumberOfPoints();
   cout << "Found "<<n<<" points in the file..." << endl;
   vtkDoubleArray *dens = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("dens"));
   vtkDoubleArray *h = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("h"));
   vtkDoubleArray *u = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("u"));
   vtkDoubleArray *mass = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("mass"));
   vtkDoubleArray *colour = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("colour"));
   vtkIntArray *iam = static_cast<vtkIntArray *> (dataset->GetPointData()->GetArray("iam"));
   vtkIntArray *tag= static_cast<vtkIntArray *> (dataset->GetPointData()->GetArray("tag"));
   vtkDoubleArray *v = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("v"));
   vtkDoubleArray *vort = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("vort"));
#ifndef BACKCOMPAT_READ
   vtkDoubleArray *vhat = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("vhat"));
   vtkDoubleArray *norm1 = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("norm1"));
   vtkDoubleArray *norm2 = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("norm2"));
#ifdef _3D_
   vtkDoubleArray *norm3 = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("norm3"));
#endif
   vtkDoubleArray *eFF = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("eFF"));
#endif
   vtkDoubleArray *alpha = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("alpha"));
   vtkDoubleArray *eBForce= static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("eBForce"));
   vtkDoubleArray *eViscF = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("eViscF"));
   vtkDoubleArray *eViscB = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("eViscB"));


   outPs->reserve(MAX_NUM_PARTICLES_PER_CPU);
   outPs->resize(n);
   
   vector<particleT>::iterator p = outPs->begin();

   for (int i=0;i<n;i++) {
      readRestartProcessTuple(i,vtkData,data,vars);
      double r[3],vv[3],nn1[3],nn2[3],vvhat[3],vvort[3];
#ifdef _3D_
      double nn3[3];
      norm3->GetTuple(i,nn3);
#endif
      dataset->GetPoint(i,r);
      v->GetTuple(i,vv);
#ifndef BACKCOMPAT_READ
      vort->GetTuple(i,vvort);
      vhat->GetTuple(i,vvhat);
      norm1->GetTuple(i,nn1);
      norm2->GetTuple(i,nn2);
#endif
      for (int j=0;j<NDIM;j++) {
         p->r[j]=r[j];
         p->v[j]=vv[j];
#ifndef BACKCOMPAT_READ
         p->vort[j]=vvort[j];
         p->vhat[j]=vvhat[j];
         p->norm1[j]=nn1[j];
         p->norm2[j]=nn2[j];
#ifdef _3D_
         p->norm3[j]=nn3[j];
#endif
#else
         p->vhat[j]=vv[j];
#endif
      }
      //cout << "Processor "<<globals->mpiRank<<": got particle at "<<p->r<<endl;
      p->dens = dens->GetValue(i);
      p->h = h->GetValue(i);
      p->eBForce = eBForce->GetValue(i);
      p->eViscF = eViscF->GetValue(i);
      p->eViscB = eViscB->GetValue(i);
#ifndef BACKCOMPAT_READ
      p->eFF = eFF->GetValue(i);
#endif
      p->u = u->GetValue(i);
      p->mass = mass->GetValue(i);
      p->colour = colour->GetValue(i);
      p->iam = (iamTypes)iam->GetValue(i);
      p->tag= tag->GetValue(i);
      p->alpha= alpha->GetValue(i);
#ifdef SLK
      p->currR = p->r;
#endif
      p++;
   }

   reader->Delete();
}

void Cio_data_vtk::readOutputDem(int timestep,vector<demParticle *outPs,CglobalVars *globals) {}

void Cio_data_vtk::readOutputSph(int timestep,vector<sphParticle *outPs,CglobalVars *globals) {
   
   vtkXMLPUnstructuredGridReader *reader = vtkXMLPUnstructuredGridReader::New();
   char strRank[20];
   sprintf(strRank,"%d",globals->mpiRank);
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   string theFilename = filename+strTimestep+".pvtu";
   reader->SetFileName(theFilename.c_str());
   reader->Update(); 
   vtkUnstructuredGrid *dataset = reader->GetOutput();
   
   int n = dataset->GetNumberOfPoints();
   //vtkFloatArray *colour = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("colour"));
   vtkIntArray *tag= static_cast<vtkIntArray *> (dataset->GetPointData()->GetArray("tag"));
   vtkIntArray *iam = static_cast<vtkIntArray *> (dataset->GetPointData()->GetArray("iam"));
   vtkFloatArray *v = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("v"));
#ifdef SAVE_VHALF
   vtkFloatArray *vHalf = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("vHalf"));
   vtkFloatArray *rHalf = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("rHalf"));
#endif
   vtkFloatArray *vhat = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("vhat"));
#ifdef AVE_VELOCITY
   vtkFloatArray *aveV = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("aveV"));
#endif
   vtkFloatArray *fv = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("fv"));
   vtkFloatArray *vort= static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("vort"));
   vtkFloatArray *h= static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("h"));
   vtkFloatArray *dens = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("dens"));
   vtkFloatArray *press = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("press"));

   outPs->clear();
   outPs->resize(n);
   
   particleContainer::iterator p = outPs->begin();

   for (int i=0;i<n;i++) {
      double r[3],vv[3],fvv[3],vvhat[3],vvort[3];
#ifdef AVE_VELOCITY
      double aaveV[3];
#endif

      dataset->GetPoint(i,r);
      v->GetTuple(i,vv);
      vhat->GetTuple(i,vvhat);
      vort->GetTuple(i,vvort);
#ifdef AVE_VELOCITY
      aveV->GetTuple(i,aaveV);
#endif
      fv->GetTuple(i,fvv);
#ifdef SAVE_VHALF
      double vvHalf[3];
      double rrHalf[3];
      vHalf->GetTuple(i,vvHalf);
      rHalf->GetTuple(i,rrHalf);
#endif
      for (int j=0;j<NDIM;j++) {
         p->r[j]=r[j];
         p->v[j]=vv[j];
#ifdef SAVE_VHALF
         p->vHalf[j]=vvHalf[j];
         p->rHalf[j]=rrHalf[j];
#endif
         p->vhat[j]=vvhat[j];
#ifdef AVE_VELOCITY
         p->aveV[j]=aaveV[j];
#endif
         p->fv[j]=fvv[j];
         p->vort[j]=vvort[j];
      }
      //p->colour = colour->GetValue(i);
      p->tag= tag->GetValue(i);
      p->iam = (iamTypes)iam->GetValue(i);
      p->h= h->GetValue(i);
      p->dens = dens->GetValue(i);
      p->press = press->GetValue(i);
      //cout <<"Reading point with r = "<<p->r<<" v = "<<p->v<<" colour = "<<p->colour<<" h = "<<p->h<<" tag = "<<p->tag<<endl;
      p++;
   }
   
   reader->Delete();
}

void Cio_data_vtk::writeRestartDem(int timestep,vector<sphParticle> &ps,CglobalVars *globals) {}
void Cio_data_vtk::writeRestartSph(int timestep,vector<demParticle> &ps,CglobalVars *globals) {
   
   //construct points

   int n = ps.size();
   vtkPoints *newPts = vtkPoints::New();;
   vtkDoubleArray *dens = vtkDoubleArray::New();
   vtkDoubleArray *h = vtkDoubleArray::New();
   vtkDoubleArray *u = vtkDoubleArray::New();
   vtkDoubleArray *mass = vtkDoubleArray::New();
   vtkDoubleArray *colour = vtkDoubleArray::New();
   vtkIntArray *iam = vtkIntArray::New();
   vtkDoubleArray *v = vtkDoubleArray::New();
   vtkDoubleArray *vhat = vtkDoubleArray::New();
   vtkDoubleArray *norm1 = vtkDoubleArray::New();
   vtkDoubleArray *norm2 = vtkDoubleArray::New();
#ifdef _3D_
   vtkDoubleArray *norm3 = vtkDoubleArray::New();
#endif
   vtkDoubleArray *alpha = vtkDoubleArray::New();
   vtkDoubleArray *vort = vtkDoubleArray::New();
   vtkIntArray *tag = vtkIntArray::New();
   vtkDoubleArray *eBForce = vtkDoubleArray::New();
   vtkDoubleArray *eViscF = vtkDoubleArray::New();
   vtkDoubleArray *eViscB = vtkDoubleArray::New();
   vtkDoubleArray *eFF = vtkDoubleArray::New();

   newPts->SetNumberOfPoints(n);

   dens->SetNumberOfValues(n);
   h->SetNumberOfValues(n);
   eBForce->SetNumberOfValues(n);
   eViscF->SetNumberOfValues(n);
   eViscB->SetNumberOfValues(n);
   eFF->SetNumberOfValues(n);
   u->SetNumberOfValues(n);
   mass->SetNumberOfValues(n);
   colour->SetNumberOfValues(n);
   iam->SetNumberOfValues(n);
   tag->SetNumberOfValues(n);
   v->SetNumberOfComponents(3);
   v->SetNumberOfTuples(n);
   vhat->SetNumberOfComponents(3);
   vhat->SetNumberOfTuples(n);
   vort->SetNumberOfComponents(3);
   vort->SetNumberOfTuples(n);
   norm1->SetNumberOfComponents(3);
   norm1->SetNumberOfTuples(n);
   norm2->SetNumberOfComponents(3);
   norm2->SetNumberOfTuples(n);
#ifdef _3D_
   norm3->SetNumberOfComponents(3);
   norm3->SetNumberOfTuples(n);
#endif
   alpha->SetNumberOfValues(n);

   dens->SetName("dens");
   h->SetName("h");
   eBForce->SetName("eBForce");
   eViscF->SetName("eViscF");
   eViscB->SetName("eViscB");
   eFF->SetName("eFF");
   u->SetName("u");
   mass->SetName("mass");
   colour->SetName("colour");
   iam->SetName("iam");
   tag->SetName("tag");
   v->SetName("v");
   vhat->SetName("vhat");
   vort->SetName("vort");
   norm1->SetName("norm1");
   norm2->SetName("norm2");
#ifdef _3D_
   norm3->SetName("norm3");
#endif
   alpha->SetName("alpha");

   int i=0;
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
#if NDIM==1
      newPts->SetPoint(i,p->r[0],0,0);
      v->SetTuple3(i,p->v[0],0,0); 
      vhat->SetTuple3(i,p->vhat[0],0,0); 
      vort->SetTuple3(i,p->vort[0],0,0);
      norm1->SetTuple3(i,p->norm1[0],0,0); 
      norm2->SetTuple3(i,p->norm2[0],0,0); 
#elif NDIM==2
      newPts->SetPoint(i,p->r[0],p->r[1],0);
      v->SetTuple3(i,p->v[0],p->v[1],0); 
      vhat->SetTuple3(i,p->vhat[0],p->vhat[1],0); 
      vort->SetTuple3(i,p->vort[0],p->vort[1],0);
      norm1->SetTuple3(i,p->norm1[0],p->norm1[1],0);
      norm2->SetTuple3(i,p->norm2[0],p->norm2[1],0);
#elif NDIM==3
      newPts->SetPoint(i,p->r[0],p->r[1],p->r[2]);
      v->SetTuple3(i,p->v[0],p->v[1],p->v[2]); 
      vhat->SetTuple3(i,p->vhat[0],p->vhat[1],p->vhat[2]); 
      vort->SetTuple3(i,p->vort[0],p->vort[1],p->vort[2]); 
      norm1->SetTuple3(i,p->norm1[0],p->norm1[1],p->norm1[2]);
      norm2->SetTuple3(i,p->norm2[0],p->norm2[1],p->norm2[2]);
      norm3->SetTuple3(i,p->norm3[0],p->norm3[1],p->norm3[2]);
#endif
      dens->SetValue(i,p->dens);
      h->SetValue(i,p->h);
      eBForce->SetValue(i,p->eBForce);
      eViscF->SetValue(i,p->eViscF);
      eViscB->SetValue(i,p->eViscB);
      eFF->SetValue(i,p->eFF);
      u->SetValue(i,p->u);
      mass->SetValue(i,p->mass);
      colour->SetValue(i,p->colour);
      iam->SetValue(i,p->iam);
      tag->SetValue(i,p->tag);
      alpha->SetValue(i,p->alpha);
      i++;
   }

   vtkUnstructuredGrid *dataset = vtkUnstructuredGrid::New();
   dataset->SetPoints(newPts);
   dataset->GetPointData()->AddArray(dens);
   dataset->GetPointData()->AddArray(h);
   dataset->GetPointData()->AddArray(eBForce);
   dataset->GetPointData()->AddArray(eViscF);
   dataset->GetPointData()->AddArray(eViscB);
   dataset->GetPointData()->AddArray(eFF);
   dataset->GetPointData()->AddArray(u);
   dataset->GetPointData()->AddArray(mass);
   dataset->GetPointData()->AddArray(colour);
   dataset->GetPointData()->AddArray(iam);
   dataset->GetPointData()->AddArray(tag);
   dataset->GetPointData()->AddArray(v);
   dataset->GetPointData()->AddArray(vhat);
   dataset->GetPointData()->AddArray(vort);
   dataset->GetPointData()->AddArray(norm1);
   dataset->GetPointData()->AddArray(norm2);
#ifdef _3D_
   dataset->GetPointData()->AddArray(norm3);
#endif
   dataset->GetPointData()->AddArray(alpha);
   dataset->GetPointData()->SetActiveScalars("iam");
   dataset->GetPointData()->SetActiveVectors("v");
   dataset->GetPointData()->SetActiveNormals("norm1");
   dataset->GetPointData()->SetActiveGlobalIds("tag");

   //write data
   vtkXMLPUnstructuredGridWriter *writer = vtkXMLPUnstructuredGridWriter::New();
   writer->SetNumberOfPieces(globals->mpiSize);
   writer->SetStartPiece(globals->mpiRank);
   writer->SetEndPiece(globals->mpiRank);
   writer->SetInput(dataset);
   writer->SetDataModeToBinary();
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   string outFilename = filename+"Restart"+strTimestep+".pvtu";
   writer->SetFileName(outFilename.c_str());
   cout << "\tWriting "<<n<<" particles to file: "<<outFilename<<endl;
   //writer->SetTimeStep(timestep);
   writer->Write();


   writer->Delete();
   dataset->Delete();
   newPts->Delete();
   dens->Delete();
   u->Delete();
   v->Delete();
   vhat->Delete();
   h->Delete();
   eBForce->Delete();
   eViscF->Delete();
   eViscB->Delete();
   eFF->Delete();
   mass->Delete();
   colour->Delete();
   iam->Delete();
   tag->Delete();
   vort->Delete();
   norm1->Delete();
   norm2->Delete();
#ifdef _3D_
   norm3->Delete();
#endif
   
   alpha->Delete();
}

void Cio_data_vtk::writeOutput(int timestep,particleContainer &ps,CcustomSimBase &customSim,CglobalVars *globals) {
   //construct points

   int n = ps.size();
   vtkPoints *newPts = vtkPoints::New();
   vtkCellArray *newCells = vtkCellArray::New();
   vtkFloatArray *dens = vtkFloatArray::New();
   vtkFloatArray *press= vtkFloatArray::New();
   //vtkFloatArray *dddt= vtkFloatArray::New();
   //vtkFloatArray *colour = vtkFloatArray::New();
   vtkFloatArray *v = vtkFloatArray::New();
   vtkFloatArray *vhat = vtkFloatArray::New();
#ifdef LIQ_DEM
   vtkFloatArray *porosity= vtkFloatArray::New();
   vtkFloatArray *fdrag = vtkFloatArray::New();
#endif

#ifdef MY_VAR_RES
   vtkFloatArray *gradVx= vtkFloatArray::New();
   vtkFloatArray *gradVy= vtkFloatArray::New();
#endif
#ifdef SLK
   vtkFloatArray *dr= vtkFloatArray::New();
   vtkFloatArray *mass= vtkFloatArray::New();
#endif
   
   
#ifdef AVE_VELOCITY
   vtkFloatArray *aveV = vtkFloatArray::New();
#endif
#ifdef SAVE_VHALF
   vtkDoubleArray *vHalf = vtkDoubleArray::New();
   vtkDoubleArray *rHalf = vtkDoubleArray::New();
#endif
   //vtkFloatArray *viscv = vtkFloatArray::New();
   //vtkFloatArray *norm1 = vtkFloatArray::New();
   //vtkFloatArray *norm2 = vtkFloatArray::New();
   vtkFloatArray *ff = vtkFloatArray::New();
   vtkFloatArray *fb = vtkFloatArray::New();
   vtkFloatArray *fp = vtkFloatArray::New();
   vtkFloatArray *fv = vtkFloatArray::New();
   vtkFloatArray *vort = vtkFloatArray::New();
   vtkFloatArray *h= vtkFloatArray::New();
   vtkIdTypeArray *tag= vtkIdTypeArray::New();
   vtkIntArray *iam = vtkIntArray::New();
   //vtkFloatArray *dudt= vtkFloatArray::New();
   //vtkFloatArray *u = vtkFloatArray::New();
   //vtkFloatArray *tmp = vtkFloatArray::New();
   //vtkFloatArray *deViscFdt = vtkFloatArray::New();
   //vtkFloatArray *deViscBdt = vtkFloatArray::New();
   //vtkFloatArray *deBForcedt = vtkFloatArray::New();

   newPts->SetNumberOfPoints(n);

   dens->SetNumberOfValues(n);
   press->SetNumberOfValues(n);
   //dddt->SetNumberOfValues(n);
   //colour->SetNumberOfValues(n);
   tag->SetNumberOfValues(n);
   iam->SetNumberOfValues(n);
#ifdef LIQ_DEM
   porosity->SetNumberOfValues(n);
   fdrag->SetNumberOfComponents(3);
   fdrag->SetNumberOfTuples(n);
#endif
#ifdef SAVE_VHALF
   vHalf->SetNumberOfComponents(3);
   vHalf->SetNumberOfTuples(n);
   rHalf->SetNumberOfComponents(3);
   rHalf->SetNumberOfTuples(n);
#endif
#ifdef SLK
   dr->SetNumberOfComponents(3);
   dr->SetNumberOfTuples(n);
   mass->SetNumberOfValues(n);
#endif
   v->SetNumberOfComponents(3);
   v->SetNumberOfTuples(n);
   vhat->SetNumberOfComponents(3);
   vhat->SetNumberOfTuples(n);
#ifdef MY_VAR_RES
   gradVx->SetNumberOfComponents(3);
   gradVx->SetNumberOfTuples(n);
   gradVy->SetNumberOfComponents(3);
   gradVy->SetNumberOfTuples(n);
#endif
#ifdef AVE_VELOCITY
   aveV->SetNumberOfComponents(3);
   aveV->SetNumberOfTuples(n);
#endif
   //viscv->SetNumberOfComponents(3);
   //viscv->SetNumberOfTuples(n);
   //norm1->SetNumberOfComponents(3);
   //norm1->SetNumberOfTuples(n);
   //norm2->SetNumberOfComponents(3);
   //norm2->SetNumberOfTuples(n);
   vort->SetNumberOfComponents(3);
   vort->SetNumberOfTuples(n);
   h->SetNumberOfValues(n);
   //dudt->SetNumberOfValues(n);
   //u->SetNumberOfValues(n);
   //tmp->SetNumberOfValues(n);
   //deViscFdt->SetNumberOfValues(n);
   //deViscBdt->SetNumberOfValues(n);
   //deBForcedt->SetNumberOfValues(n);
   ff->SetNumberOfComponents(3);
   ff->SetNumberOfTuples(n);
   fb->SetNumberOfComponents(3);
   fb->SetNumberOfTuples(n);
   fv->SetNumberOfComponents(3);
   fv->SetNumberOfTuples(n);
   fp->SetNumberOfComponents(3);
   fp->SetNumberOfTuples(n);

   dens->SetName("dens");
   press->SetName("press");
#ifdef LIQ_DEM
   porosity->SetName("porosity");
   fdrag->SetName("fdrag");
#endif
   //dddt->SetName("dddt");
   //colour->SetName("colour");
#ifdef SAVE_VHALF
   vHalf->SetName("vHalf");
   rHalf->SetName("rHalf");
#endif
#ifdef SLK
   dr->SetName("dr");
   mass->SetName("mass");
#endif

   v->SetName("v");
   vhat->SetName("vhat");
#ifdef MY_VAR_RES
   gradVx->SetName("gradVx");
   gradVy->SetName("gradVy");
#endif
#ifdef AVE_VELOCITY
   aveV->SetName("aveV");
#endif
   //viscv->SetName("viscv");
   //norm1->SetName("norm1");
   //norm2->SetName("norm2");
   vort->SetName("vort");
   h->SetName("h");
   tag->SetName("tag");
   iam->SetName("iam");
   //dudt->SetName("dudt");
   //u->SetName("u");
   //tmp->SetName("tmp");
   //deViscFdt->SetName("deViscFdt");
   //deViscBdt->SetName("deViscBdt");
   //deBForcedt->SetName("deBForcedt");
   ff->SetName("ff");
   fb->SetName("fb");
   fv->SetName("fv");
   fp->SetName("fp");

   int i = 0;
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
      newCells->InsertNextCell(1);
      newCells->InsertCellPoint(i);
      vect newr = customSim.rFilter(p->r);
      vect newv = customSim.vFilter(p->v,p->r);
      vect newvhat = customSim.vFilter(p->vhat,p->r);
#ifdef MY_VAR_RES
      vect newgradVx = customSim.vFilter(p->gradV[0],p->r);
      vect newgradVy = customSim.vFilter(p->gradV[1],p->r);
#endif
#ifdef AVE_VELOCITY
      vect newaveV = customSim.vFilter(p->aveV,p->r);
#endif
#ifdef SAVE_VHALF
      vect newrHalf = customSim.rFilter(p->rHalf);
      vect newvHalf = customSim.vFilter(p->vHalf,p->rHalf);
#endif
      //vect newviscv = customSim.vFilter(p->viscv,p->r);
#if NDIM==1
      newPts->SetPoint(i,newr[0],0,0);
      v->SetTuple3(i,newv[0],0,0); 
      vhat->SetTuple3(i,newvhat[0],0,0); 
#ifdef AVE_VELOCITY
      aveV->SetTuple3(i,newaveV[0],0,0); 
#endif
      //norm1->SetTuple3(i,p->norm1[0],0,0); 
      //norm2->SetTuple3(i,p->norm2[0],0,0); 
      vort->SetTuple3(i,p->vort[0],0,0);
      ff->SetTuple3(i,p->ff[0],0,0); 
      fb->SetTuple3(i,p->fb[0],0,0); 
      fv->SetTuple3(i,p->fv[0],0,0); 
      fp->SetTuple3(i,p->fp[0],0,0); 
#elif NDIM==2
      newPts->SetPoint(i,newr[0],newr[1],0);
#ifdef SAVE_VHALF
      vHalf->SetTuple3(i,newvHalf[0],newvHalf[1],0); 
      rHalf->SetTuple3(i,newrHalf[0],newrHalf[1],0); 
#endif
      v->SetTuple3(i,newv[0],newv[1],0); 
      vhat->SetTuple3(i,newvhat[0],newvhat[1],0); 
#ifdef MY_VAR_RES
      gradVx->SetTuple3(i,newgradVx[0],newgradVx[1],0); 
      gradVy->SetTuple3(i,newgradVy[0],newgradVy[1],0); 
#endif
#ifdef SLK
      vect tmp = p->currR-p->r;
      dr->SetTuple3(i,tmp[0],tmp[1],0);
#endif
#ifdef AVE_VELOCITY
      aveV->SetTuple3(i,newaveV[0],newaveV[1],0); 
#endif
      vort->SetTuple3(i,p->vort[0],p->vort[1],0);
      ff->SetTuple3(i,p->ff[0],p->ff[1],0); 
      fb->SetTuple3(i,p->fb[0],p->fb[1],0); 
      fv->SetTuple3(i,p->fv[0],p->fv[1],0); 
      fp->SetTuple3(i,p->fp[0],p->fp[1],0); 
#elif NDIM==3
      newPts->SetPoint(i,newr[0],newr[1],newr[2]);
      v->SetTuple3(i,newv[0],newv[1],newv[2]); 
      vhat->SetTuple3(i,newvhat[0],newvhat[1],newvhat[2]); 
#ifdef AVE_VELOCITY
      aveV->SetTuple3(i,newaveV[0],newaveV[1],newaveV[2]); 
#endif
      vort->SetTuple3(i,p->vort[0],p->vort[1],p->vort[2]); 
      ff->SetTuple3(i,p->ff[0],p->ff[1],p->ff[2]); 
      fb->SetTuple3(i,p->fb[0],p->fb[1],p->fb[2]); 
      fv->SetTuple3(i,p->fv[0],p->fv[1],p->fv[2]); 
      fp->SetTuple3(i,p->fp[0],p->fp[1],p->fp[2]); 
#endif
#ifdef SLK
      mass->SetValue(i,p->mass);
#endif
#ifdef LIQ_DEM
      porosity->SetValue(i,p->porosity);
#if NDIM==2
      fdrag->SetTuple3(i,p->fdrag[0],p->fdrag[1],0.0);
#else
      fdrag->SetTuple3(i,p->fdrag[0],p->fdrag[1],p->fdrag[2]);
#endif

#endif
      dens->SetValue(i,p->dens);
      press->SetValue(i,p->press);
      //dddt->SetValue(i,p->dddt);
      //colour->SetValue(i,p->colour);
      h->SetValue(i,p->h);
      tag->SetValue(i,p->tag);
      iam->SetValue(i,p->iam);
      i++;
   }

   vtkUnstructuredGrid *dataset = vtkUnstructuredGrid::New();
   dataset->SetPoints(newPts);
   dataset->SetCells(1,newCells);
   dataset->GetPointData()->AddArray(dens);
   dataset->GetPointData()->AddArray(press);
   //dataset->GetPointData()->AddArray(dddt);
   //dataset->GetPointData()->AddArray(colour);
   dataset->GetPointData()->AddArray(tag);
   dataset->GetPointData()->AddArray(iam);
#ifdef SAVE_VHALF
   dataset->GetPointData()->AddArray(vHalf);
   dataset->GetPointData()->AddArray(rHalf);
#endif
   dataset->GetPointData()->AddArray(v);
   dataset->GetPointData()->AddArray(vhat);
#ifdef MY_VAR_RES
   dataset->GetPointData()->AddArray(gradVx);
   dataset->GetPointData()->AddArray(gradVy);
#endif
#ifdef LIQ_DEM
   dataset->GetPointData()->AddArray(porosity);
   dataset->GetPointData()->AddArray(fdrag);
#endif
#ifdef AVE_VELOCITY
   dataset->GetPointData()->AddArray(aveV);
#endif
#ifdef SLK
   dataset->GetPointData()->AddArray(dr);
   dataset->GetPointData()->AddArray(mass);
#endif
   //dataset->GetPointData()->AddArray(viscv);
   //dataset->GetPointData()->AddArray(norm1);
   //dataset->GetPointData()->AddArray(norm2);
   dataset->GetPointData()->AddArray(vort);
   dataset->GetPointData()->AddArray(h);
   dataset->GetPointData()->AddArray(ff);
   dataset->GetPointData()->AddArray(fb);
   dataset->GetPointData()->AddArray(fv);
   dataset->GetPointData()->AddArray(fp);
   //dataset->GetPointData()->AddArray(dudt);
   //dataset->GetPointData()->AddArray(u);
   //dataset->GetPointData()->AddArray(tmp);
   //dataset->GetPointData()->AddArray(deViscFdt);
   //dataset->GetPointData()->AddArray(deViscBdt);
   //dataset->GetPointData()->AddArray(deBForcedt);
   dataset->GetPointData()->SetActiveScalars("dens");
   dataset->GetPointData()->SetActiveVectors("v");
   dataset->GetPointData()->SetActiveGlobalIds("tag");

   //write data
   vtkXMLPUnstructuredGridWriter *writer = vtkXMLPUnstructuredGridWriter::New();
   writer->SetNumberOfPieces(globals->mpiSize);
   writer->SetStartPiece(globals->mpiRank);
   writer->SetEndPiece(globals->mpiRank);
   writer->SetNumberOfTimeSteps(OUTSTEP);
   writer->SetTimeStepRange(globals->outstep,globals->outstep);
   writer->SetInput(dataset);
   writer->SetDataModeToBinary();
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   string outFilename = filename+strTimestep+".pvtu";
   writer->SetFileName(outFilename.c_str());
   if (globals->mpiRank==0) cout << "\tWriting "<<n<<" particles to file: "<<outFilename<<endl;
   writer->SetTimeStep(timestep);
   writer->WriteNextTime(globals->time);
   //writer->Write();


   writer->Delete();
   dataset->Delete();
   newPts->Delete();
   newCells->Delete();
   dens->Delete();
   press->Delete();
   //dddt->Delete();
   //dudt->Delete();
   //u->Delete();
   //tmp->Delete();
   //deViscFdt->Delete();
   //deViscBdt->Delete();
   //deBForcedt->Delete();
#ifdef LIQ_DEM
   porosity->Delete();
   fdrag->Delete();
#endif
#ifdef SAVE_VHALF
   vHalf->Delete();
   rHalf->Delete();
#endif
   v->Delete();
   vhat->Delete();
#ifdef MY_VAR_RES
   gradVx->Delete();
   gradVy->Delete();
#endif
#ifdef AVE_VELOCITY
   aveV->Delete();
#endif
   //viscv->Delete();
   //norm1->Delete();
   //norm2->Delete();
   //colour->Delete();
   vort->Delete();
   h->Delete();
   ff->Delete();
   fb->Delete();
   fv->Delete();
   fp->Delete();
   tag->Delete();
   iam->Delete();
}

void Cio_data_vtk::writeAuxNoData(int timestep,particleContainer &ps,const char *name,CglobalVars *globals) {
   int n = ps.size();
   vtkPoints *newPts = vtkPoints::New();;
   vtkIdTypeArray *tag= vtkIdTypeArray::New();

   newPts->SetNumberOfPoints(n);

   tag->SetNumberOfValues(n);

   tag->SetName("tag");

   int i = 0;
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
      vect newr = p->r;
#if NDIM==1
      newPts->SetPoint(i,newr[0],0,0);
#elif NDIM==2
      newPts->SetPoint(i,newr[0],newr[1],0);
#elif NDIM==3
      newPts->SetPoint(i,newr[0],newr[1],newr[2]);
#endif
      tag->SetValue(i,p->tag);
      i++;
   }

   vtkUnstructuredGrid *dataset = vtkUnstructuredGrid::New();
   dataset->SetPoints(newPts);
   dataset->GetPointData()->AddArray(tag);
   dataset->GetPointData()->SetActiveScalars(name);
   //dataset->GetPointData()->SetActiveGlobalIds("tag");

   //write data
   vtkXMLPUnstructuredGridWriter *writer = vtkXMLPUnstructuredGridWriter::New();
   writer->SetNumberOfPieces(globals->mpiSize);
   writer->SetStartPiece(globals->mpiRank);
   writer->SetEndPiece(globals->mpiRank);
   writer->SetInput(dataset);
   writer->SetDataModeToBinary();
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   string outFilename = filename+name+strTimestep+".pvtu";
   writer->SetFileName(outFilename.c_str());
   if (globals->mpiRank==0) cout << "\tWriting "<<n<<" particles to file: "<<outFilename<<endl;
   //writer->SetTimeStep(timestep);
   writer->Write();


   writer->Delete();
   dataset->Delete();
   newPts->Delete();
   tag->Delete();
}

void Cio_data_vtk::writeAux(int timestep,particleContainer &ps,vector<double> theData,const char *name,CglobalVars *globals) {
   int n = ps.size();
   vtkPoints *newPts = vtkPoints::New();;
   vtkFloatArray *data = vtkFloatArray::New();
   vtkIdTypeArray *tag= vtkIdTypeArray::New();

   newPts->SetNumberOfPoints(n);

   data->SetNumberOfValues(n);
   tag->SetNumberOfValues(n);

   data->SetName(name);
   tag->SetName("tag");

   int i = 0;
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
      vect newr = p->r;
#if NDIM==1
      newPts->SetPoint(i,newr[0],0,0);
#elif NDIM==2
      newPts->SetPoint(i,newr[0],newr[1],0);
#elif NDIM==3
      newPts->SetPoint(i,newr[0],newr[1],newr[2]);
#endif
      data->SetValue(i,theData[i]);
      tag->SetValue(i,p->tag);
      i++;
   }

   vtkUnstructuredGrid *dataset = vtkUnstructuredGrid::New();
   dataset->SetPoints(newPts);
   dataset->GetPointData()->AddArray(data);
   dataset->GetPointData()->AddArray(tag);
   dataset->GetPointData()->SetActiveScalars(name);
   //dataset->GetPointData()->SetActiveGlobalIds("tag");

   //write data
   vtkXMLPUnstructuredGridWriter *writer = vtkXMLPUnstructuredGridWriter::New();
   writer->SetNumberOfPieces(globals->mpiSize);
   writer->SetStartPiece(globals->mpiRank);
   writer->SetEndPiece(globals->mpiRank);
   writer->SetInput(dataset);
   writer->SetDataModeToBinary();
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   string outFilename = filename+name+strTimestep+".pvtu";
   writer->SetFileName(outFilename.c_str());
   if (globals->mpiRank==0) cout << "\tWriting "<<n<<" particles to file: "<<outFilename<<endl;
   //writer->SetTimeStep(timestep);
   writer->Write();


   writer->Delete();
   dataset->Delete();
   newPts->Delete();
   data->Delete();
   tag->Delete();
}


void Cio_data_vtk::writeGrid(int timestep,particleContainer &ps,vectInt &gridDims,CglobalVars *globals) {
   //construct points

   int n = ps.size();
   vtkPoints *newPts = vtkPoints::New();;
   vtkDoubleArray *dens = vtkDoubleArray::New();
   //vtkDoubleArray *u = vtkDoubleArray::New();
   vtkDoubleArray *v = vtkDoubleArray::New();
   vtkDoubleArray *vort = vtkDoubleArray::New();

   newPts->SetNumberOfPoints(n);

   dens->SetNumberOfValues(n);
   //u->SetNumberOfValues(n);
   v->SetNumberOfComponents(3);
   v->SetNumberOfTuples(n);
   if (NDIM==3) {
      vort->SetNumberOfComponents(3);
      vort->SetNumberOfTuples(n);
   } else {
      vort->SetNumberOfValues(n);
   }

   dens->SetName("dens");
   //u->SetName("u");
   v->SetName("v");
   vort->SetName("vort");

   int i=0; 
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
#if NDIM==1
      newPts->SetPoint(i,p->r[0],0,0);
      vort->SetValue(i,p->vort);
      v->SetTuple3(i,p->v[0],0,0);
#elif NDIM==2 
      newPts->SetPoint(i,p->r[0],p->r[1],0);
      vort->SetValue(i,p->vort[0]);
      v->SetTuple3(i,p->v[0],p->v[1],0);
#elif NDIM==3 
      newPts->SetPoint(i,p->r[0],p->r[1],p->r[2]);
      v->SetTuple3(i,p->v[0],p->v[1],p->v[2]);
      vort->SetTuple3(i,p->vort[0],p->vort[1],p->vort[2]);
#endif
      dens->SetValue(i,p->dens);
      //u->SetValue(i,p->u);
      i++;
   }

   vtkStructuredGrid *dataset = vtkStructuredGrid::New();
   switch (NDIM) {
      case 1: dataset->SetDimensions(gridDims[0],1,1); 
              dataset->SetExtent(0,gridDims[0]-1,0,0,0,0); 
              break;
      case 2: dataset->SetDimensions(gridDims[0],gridDims[1],1); 
              dataset->SetExtent(0,gridDims[0]-1,0,gridDims[1]-1,0,0); 
              break;
      case 3: dataset->SetDimensions(gridDims[0],gridDims[1],gridDims[2]); 
              dataset->SetExtent(0,gridDims[0]-1,0,gridDims[1]-1,0,gridDims[2]-1); 
              break;
   }
   dataset->SetPoints(newPts);
   dataset->GetPointData()->AddArray(dens);
   //dataset->GetPointData()->AddArray(u);
   dataset->GetPointData()->AddArray(v);
   dataset->GetPointData()->AddArray(vort);
   dataset->GetPointData()->SetActiveScalars("dens");
   dataset->GetPointData()->SetActiveVectors("v");

   //write data
   vtkXMLStructuredGridWriter *writer = vtkXMLStructuredGridWriter::New();
   writer->SetInput(dataset);
   writer->SetDataModeToBinary();
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   string outFilename = filename+"Grid"+strTimestep+".pvts";
   writer->SetFileName(outFilename.c_str());
   if (globals->mpiRank==0) cout << "\tWriting "<<n<<" grid points to file: "<<outFilename<<endl;
   //writer->SetTimeStep(timestep);
   writer->Write();

   writer->Delete();
   dataset->Delete();
   newPts->Delete();
   dens->Delete();
   //u->Delete();
   v->Delete();
   vort->Delete();
}
