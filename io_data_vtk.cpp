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

void Cio_data_vtk::readRestart(int timestep,particleContainer *outPs,CglobalVars *globals) {
   vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();

   //setup filename from given base filename
   char strTimestep[20];
   sprintf(strTimestep,"%7.7d",timestep);
   char strRank[20];
   sprintf(strRank,"%d",globals->mpiRank);
   string theFilename = filename+"Restart"+strTimestep+"_"+strRank+".vtu";
   cout << "Reading file " << theFilename << endl;
   reader->SetFileName(theFilename.c_str());

   //read data file
   reader->Update(); 
   vtkUnstructuredGrid *dataset = reader->GetOutput();
   
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
   vtkDoubleArray *eBForce= static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("eBForce"));
   vtkDoubleArray *eViscF = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("eViscF"));
   vtkDoubleArray *eViscB = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("eViscB"));

   outPs->reserve(MAX_NUM_PARTICLES_PER_CPU);
   outPs->resize(n);
   
   particleContainer::iterator p = outPs->begin();

   for (int i=0;i<n;i++) {
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
      //p->alpha= alpha->GetValue(i);
#ifdef SLK
      p->currR = p->r;
#endif
      p++;
   }

   reader->Delete();
}

void Cio_data_vtk::readOutput(int timestep,particleContainer *outPs,CglobalVars *globals) {
   
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
   vtkIntArray *tag= static_cast<vtkIntArray *> (dataset->GetPointData()->GetArray("tag"));
   vtkIntArray *iam = static_cast<vtkIntArray *> (dataset->GetPointData()->GetArray("iam"));
#ifdef LIQ_DEM
   vtkFloatArray *shepSum = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("shepSum"));
   vtkFloatArray *porosity = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("porosity"));
#endif
   vtkFloatArray *v = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("v"));
   vtkFloatArray *vhat = static_cast<vtkFloatArray *> (dataset->GetPointData()->GetArray("vhat"));
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

      dataset->GetPoint(i,r);
      v->GetTuple(i,vv);
      vhat->GetTuple(i,vvhat);
      vort->GetTuple(i,vvort);
      fv->GetTuple(i,fvv);
      for (int j=0;j<NDIM;j++) {
         p->r[j]=r[j];
         p->v[j]=vv[j];
         p->vhat[j]=vvhat[j];
         p->fv[j]=fvv[j];
         p->vort[j]=vvort[j];
      }
      p->tag= tag->GetValue(i);
#ifdef LIQ_DEM
      p->shepSum = shepSum->GetValue(i);
      p->porosity = porosity->GetValue(i);
#endif
      p->mass = pow(PSEP,NDIM)*DENS;
      p->iam = (iamTypes)iam->GetValue(i);
      p->h= h->GetValue(i);
      p->dens = dens->GetValue(i);
      p->press = press->GetValue(i);
      p++;
   }
   
   reader->Delete();
}

void Cio_data_vtk::readCSIRO(int timestep,particleContainer *outPs,CglobalVars *globals) {
   
   vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
   char strRank[20];
   sprintf(strRank,"%d",globals->mpiRank);
   char strTimestep[20];
   sprintf(strTimestep,"%5.5d",timestep);
   string theFilename = filename+strTimestep+".vtk";
   reader->SetFileName(theFilename.c_str());
   vtkUnstructuredGrid *dataset = vtkUnstructuredGrid::New();
   reader->SetOutput(dataset);
   reader->ReadAllScalarsOn();
   reader->Update(); 

   int n = dataset->GetNumberOfPoints();
   cout <<"dataset contains "<<n<<" points..."<<endl;
   vtkLongArray *tag= static_cast<vtkLongArray *> (dataset->GetPointData()->GetArray("ID"));
   //vtkDoubleArray *dens = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("density"));
   vtkDoubleArray *vx = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("vx"));
   vtkDoubleArray *vy = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("vy"));
   vtkDoubleArray *vz = static_cast<vtkDoubleArray *> (dataset->GetPointData()->GetArray("vz"));

   outPs->clear();
   outPs->resize(n);
   
   particleContainer::iterator p = outPs->begin();

   for (int i=0;i<n;i++) {
      double r[3];
      dataset->GetPoint(i,r);
      for (int j=0;j<NDIM;j++) {
         p->r[j]=r[j];
      }
      p->tag= tag->GetValue(i);
      //p->dens = dens->GetValue(i);
      p->dens = DENS;
      p->v[0]=vx->GetValue(i);
      p->v[1]=vy->GetValue(i);
      p->v[2]=vz->GetValue(i);
      p->vhat=p->v;
      p->iam = sph;
      p->mass = pow(PSEP,NDIM)*DENS;
      p->h = HFAC*PSEP;
      p++;
   }
   
   reader->Delete();
}



void Cio_data_vtk::writeRestart(int timestep,particleContainer &ps,CglobalVars *globals) {
   int n = ps.size();
   vtkCellArray *newCells = vtkCellArray::New();
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

   int i=0;
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
      newCells->InsertNextCell(1);
      newCells->InsertCellPoint(i);
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
      i++;
   }

   vtkUnstructuredGrid *dataset = vtkUnstructuredGrid::New();
   dataset->SetPoints(newPts);
   dataset->SetCells(1,newCells);
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
   //dataset->GetPointData()->AddArray(alpha);
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
   newCells->Delete();
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
}

void Cio_data_vtk::writeOutput(int timestep,particleContainer &ps,CglobalVars *globals) {
   int n = ps.size();
   vtkPoints *newPts = vtkPoints::New();
   vtkCellArray *newCells = vtkCellArray::New();
   vtkFloatArray *dens = vtkFloatArray::New();
   vtkFloatArray *press= vtkFloatArray::New();
   vtkFloatArray *v = vtkFloatArray::New();
   vtkFloatArray *vhat = vtkFloatArray::New();
#ifdef LIQ_DEM
   vtkFloatArray *porosity= vtkFloatArray::New();
   vtkFloatArray *dporositydt= vtkFloatArray::New();
   vtkFloatArray *fdrag = vtkFloatArray::New();
   vtkFloatArray *shepSum = vtkFloatArray::New();
#endif
#ifdef HALLER_LCS
   vtkFloatArray *eval1 = vtkFloatArray::New();
   vtkFloatArray *eval2 = vtkFloatArray::New();
   vtkFloatArray *evec1 = vtkFloatArray::New();
   vtkFloatArray *Z_surface = vtkFloatArray::New();
#endif
#ifdef VAR_H_CORRECTION2
   vtkFloatArray *gradH = vtkFloatArray::New();
#endif


#ifdef SLK
   vtkFloatArray *dr= vtkFloatArray::New();
   vtkFloatArray *mass= vtkFloatArray::New();
#endif
   
#ifdef TEST_VISC
   vtkFloatArray *eViscF= vtkFloatArray::New();
#endif

   vtkFloatArray *ff = vtkFloatArray::New();
   vtkFloatArray *fb = vtkFloatArray::New();
   vtkFloatArray *fp = vtkFloatArray::New();
   vtkFloatArray *fv = vtkFloatArray::New();
   vtkFloatArray *vort = vtkFloatArray::New();
   vtkFloatArray *h= vtkFloatArray::New();
   vtkIdTypeArray *tag= vtkIdTypeArray::New();
   vtkIntArray *iam = vtkIntArray::New();

   newPts->SetNumberOfPoints(n);

   dens->SetNumberOfValues(n);
   press->SetNumberOfValues(n);
   tag->SetNumberOfValues(n);
   iam->SetNumberOfValues(n);
#ifdef LIQ_DEM
   porosity->SetNumberOfValues(n);
   dporositydt->SetNumberOfValues(n);
   shepSum->SetNumberOfValues(n);
   fdrag->SetNumberOfComponents(3);
   fdrag->SetNumberOfTuples(n);
#endif
#ifdef HALLER_LCS
   eval1->SetNumberOfValues(n);
   eval2->SetNumberOfValues(n);
   evec1->SetNumberOfComponents(3);
   evec1->SetNumberOfTuples(n);
   Z_surface->SetNumberOfValues(n);
#endif
#ifdef SLK
   dr->SetNumberOfComponents(3);
   dr->SetNumberOfTuples(n);
   mass->SetNumberOfValues(n);
#endif
#ifdef VAR_H_CORRECTION2
   gradH->SetNumberOfComponents(3);
   gradH->SetNumberOfTuples(n);
#endif
   v->SetNumberOfComponents(3);
   v->SetNumberOfTuples(n);
   vhat->SetNumberOfComponents(3);
   vhat->SetNumberOfTuples(n);
   vort->SetNumberOfComponents(3);
   vort->SetNumberOfTuples(n);
   h->SetNumberOfValues(n);
#ifdef TEST_VISC
   eViscF->SetNumberOfValues(n);
#endif

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
   dporositydt->SetName("dporositydt");
   shepSum->SetName("shepSum");
   fdrag->SetName("fdrag");
#endif
#ifdef HALLER_LCS
   eval1->SetName("eval1");
   eval2->SetName("eval2");
   evec1->SetName("evec1");
   Z_surface->SetName("Z_surface");
#endif

#ifdef VAR_H_CORRECTION2
   gradH->SetName("gradH");
#endif
#ifdef SLK
   dr->SetName("dr");
   mass->SetName("mass");
#endif

   v->SetName("v");
   vhat->SetName("vhat");
   vort->SetName("vort");
   h->SetName("h");
   tag->SetName("tag");
   iam->SetName("iam");

#ifdef TEST_VISC
   eViscF->SetName("eViscF");
#endif
   ff->SetName("ff");
   fb->SetName("fb");
   fv->SetName("fv");
   fp->SetName("fp");

   int i = 0;
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
      newCells->InsertNextCell(1);
      newCells->InsertCellPoint(i);
#ifdef TEST_VISC
      vect newr = p->origPos;
#else
      vect newr = p->r;
#endif
      vect newv = p->v;
      vect newvhat = p->vhat;
#if NDIM==1
      newPts->SetPoint(i,newr[0],0,0);
      v->SetTuple3(i,newv[0],0,0); 
      vhat->SetTuple3(i,newvhat[0],0,0); 
      vort->SetTuple3(i,p->vort[0],0,0);
      ff->SetTuple3(i,p->ff[0],0,0); 
      fb->SetTuple3(i,p->fb[0],0,0); 
      fv->SetTuple3(i,p->fv[0],0,0); 
      fp->SetTuple3(i,p->fp[0],0,0); 
#elif NDIM==2
      newPts->SetPoint(i,newr[0],newr[1],0);
      v->SetTuple3(i,newv[0],newv[1],0); 
      vhat->SetTuple3(i,newvhat[0],newvhat[1],0); 
#ifdef SLK
      vect tmp = p->currR-p->r;
      dr->SetTuple3(i,tmp[0],tmp[1],0);
#endif
      vort->SetTuple3(i,p->vort[0],p->vort[1],0);
      ff->SetTuple3(i,p->ff[0],p->ff[1],0); 
      fb->SetTuple3(i,p->fb[0],p->fb[1],0); 
      fv->SetTuple3(i,p->fv[0],p->fv[1],0); 
      fp->SetTuple3(i,p->fp[0],p->fp[1],0); 
#ifdef VAR_H_CORRECTION2
      gradH->SetTuple3(i,p->gradH[0],p->gradH[1],0);
#endif
#elif NDIM==3
      newPts->SetPoint(i,newr[0],newr[1],newr[2]);
      v->SetTuple3(i,newv[0],newv[1],newv[2]); 
      vhat->SetTuple3(i,newvhat[0],newvhat[1],newvhat[2]); 
      vort->SetTuple3(i,p->vort[0],p->vort[1],p->vort[2]); 
      ff->SetTuple3(i,p->ff[0],p->ff[1],p->ff[2]); 
      fb->SetTuple3(i,p->fb[0],p->fb[1],p->fb[2]); 
      fv->SetTuple3(i,p->fv[0],p->fv[1],p->fv[2]); 
      fp->SetTuple3(i,p->fp[0],p->fp[1],p->fp[2]); 
#ifdef VAR_H_CORRECTION2
      gradH->SetTuple3(i,p->gradH[0],p->gradH[1],p->gradH[2]);
#endif
#endif
#ifdef SLK
      mass->SetValue(i,p->mass);
#endif
#ifdef LIQ_DEM
      porosity->SetValue(i,p->porosity);
      dporositydt->SetValue(i,p->dporositydt);
      shepSum->SetValue(i,p->shepSum);
#if NDIM==2
      fdrag->SetTuple3(i,p->fdrag[0],p->fdrag[1],0.0);
#else
      fdrag->SetTuple3(i,p->fdrag[0],p->fdrag[1],p->fdrag[2]);
#endif

#endif
#ifdef HALLER_LCS
      eval1->SetValue(i,p->eval1);
      eval2->SetValue(i,p->eval2);
#if NDIM==2
      evec1->SetTuple3(i,p->evec1[0],p->evec1[1],0); 
#elif NDIM==3
      evec1->SetTuple3(i,p->evec1[0],p->evec1[1],p->evec1[2]); 
#endif
      Z_surface->SetValue(i,p->Z_surface);
#endif
#ifdef TEST_VISC
      eViscF->SetValue(i,p->eViscF);
#endif
      dens->SetValue(i,p->dens);
      press->SetValue(i,p->press);
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
   dataset->GetPointData()->AddArray(tag);
   dataset->GetPointData()->AddArray(iam);
   dataset->GetPointData()->AddArray(v);
   dataset->GetPointData()->AddArray(vhat);
#ifdef LIQ_DEM
   dataset->GetPointData()->AddArray(porosity);
   dataset->GetPointData()->AddArray(dporositydt);
   dataset->GetPointData()->AddArray(fdrag);
   dataset->GetPointData()->AddArray(shepSum);
#endif
#ifdef HALLER_LCS
   dataset->GetPointData()->AddArray(eval1);
   dataset->GetPointData()->AddArray(eval2);
   dataset->GetPointData()->AddArray(evec1);
   dataset->GetPointData()->AddArray(Z_surface);
#endif
#ifdef VAR_H_CORRECTION2
   dataset->GetPointData()->AddArray(gradH);
#endif
#ifdef SLK
   dataset->GetPointData()->AddArray(dr);
   dataset->GetPointData()->AddArray(mass);
#endif
   dataset->GetPointData()->AddArray(vort);
   dataset->GetPointData()->AddArray(h);
#ifdef TEST_VISC
   dataset->GetPointData()->AddArray(eViscF);
#endif
   dataset->GetPointData()->AddArray(ff);
   dataset->GetPointData()->AddArray(fb);
   dataset->GetPointData()->AddArray(fv);
   dataset->GetPointData()->AddArray(fp);
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


   writer->Delete();
   dataset->Delete();
   newPts->Delete();
   newCells->Delete();
   dens->Delete();
   press->Delete();
#ifdef LIQ_DEM
   porosity->Delete();
   dporositydt->Delete();
   shepSum->Delete();
   fdrag->Delete();
#endif
#ifdef HALLER_LCS
   eval1->Delete();
   eval2->Delete();
   evec1->Delete();
   Z_surface->Delete();
#endif
#ifdef VAR_H_CORRECTION2
   gradH->Delete();
#endif
   v->Delete();
   vhat->Delete();
   vort->Delete();
   h->Delete();
#ifdef TEST_VISC
   eViscF->Delete();
#endif
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
   vtkCellArray *newCells = vtkCellArray::New();
   vtkFloatArray *data = vtkFloatArray::New();
   vtkIdTypeArray *tag= vtkIdTypeArray::New();

   newPts->SetNumberOfPoints(n);

   data->SetNumberOfValues(n);
   tag->SetNumberOfValues(n);

   data->SetName(name);
   tag->SetName("tag");

   int i = 0;
   for (particleContainer::iterator p = ps.begin();p!=ps.end();p++) {
      newCells->InsertNextCell(1);
      newCells->InsertCellPoint(i);
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
   dataset->SetCells(1,newCells);
   dataset->GetPointData()->AddArray(data);
   dataset->GetPointData()->AddArray(tag);

   dataset->GetPointData()->SetActiveScalars(name);
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
   string outFilename = filename+name+strTimestep+".pvtu";
   writer->SetFileName(outFilename.c_str());
   if (globals->mpiRank==0) cout << "\tWriting "<<n<<" particles to file: "<<outFilename<<endl;
   writer->SetTimeStep(timestep);
   writer->WriteNextTime(globals->time);




   writer->Delete();
   newCells->Delete();
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
   vtkDoubleArray *v = vtkDoubleArray::New();
#ifdef LIQ_DEM
   vtkDoubleArray *porosity = vtkDoubleArray::New();
   vtkDoubleArray *shepSum = vtkDoubleArray::New();
#endif

   vtkDoubleArray *vort = vtkDoubleArray::New();

   newPts->SetNumberOfPoints(n);

   dens->SetNumberOfValues(n);
#ifdef LIQ_DEM
   porosity->SetNumberOfValues(n);
   shepSum->SetNumberOfValues(n);
#endif
   v->SetNumberOfComponents(3);
   v->SetNumberOfTuples(n);
   if (NDIM==3) {
      vort->SetNumberOfComponents(3);
      vort->SetNumberOfTuples(n);
   } else {
      vort->SetNumberOfValues(n);
   }

   dens->SetName("dens");
#ifdef LIQ_DEM
   porosity->SetName("porosity");
   shepSum->SetName("shepSum");
#endif
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
#ifdef LIQ_DEM
      porosity->SetValue(i,p->porosity);
      shepSum->SetValue(i,p->shepSum);
#endif
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
#ifdef LIQ_DEM
   dataset->GetPointData()->AddArray(porosity);
   dataset->GetPointData()->AddArray(shepSum);
#endif
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
#ifdef LIQ_DEM
   porosity->Delete();
   shepSum->Delete();
#endif
   v->Delete();
   vort->Delete();
}
