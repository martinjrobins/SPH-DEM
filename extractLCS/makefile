CPPC = mpicxx
CODEDIR = ..
INCL = -I$(HOME)/usr/bin    -I. -I/opt/sw/gsl-1.10/include -I/opt/sw/vtk-5.0.3/include/vtk-5.0 -I/usr/include/vtk-5.2 -I/usr/include/vtk-5.8 -I/usr/include/vtk-5.4 -I$(HOME)/usr/include/vtk-5.8 -I$(HOME)/usr/include -I$(CODEDIR)
CPPFLAGS =  -Wno-deprecated $(INCL) -O3
#CPPFLAGS =  -Wno-deprecated $(INCL) -g
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkDICOMParser -lvtkpng -lvtkjpeg -lvtktiff -lvtksys -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
#VTKLIBS =  -lexpat -lxml2  -lvtkzlib -lvtklibxml2 -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -ldl -lvtksys
VTKLIBS =  -lexpat  -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -ldl -lvtksys
FFTWLIBS = -lfftw3 -lm
GSLLIBS = -lgsl -lgslcblas
LIBS = -L$(HOME)/usr/lib -L$(HOME)/usr/lib/vtk-5.8 -L/usr/mpi/gcc/openmpi-1.4.3/lib64/ -L/opt/sw/gsl-1.10/lib -L/opt/sw/vtk-5.0.3/lib $(LDFLAGS) $(VTKLIBS) $(GSLLIBS)


%.o : $(CODEDIR)/%.cpp
	$(CPPC) $(CPPFLAGS) -c $< -o $@

%.o : %.cpp
	$(CPPC) $(CPPFLAGS) -c $< -o $@

NONSPH = io_data_vtk.o dataLL.o vect.o ioGlobals.o ioDomain.o misc.o

extractLCS: $(NONSPH) extractLCS.o
	$(CPPC) $(CPPFLAGS) -o extractLCS $(NONSPH) extractLCS.o $(LIBS)
clean:
	rm *.o extractLCS 
