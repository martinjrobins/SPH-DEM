CPPC = mpicxx
CODEDIR = ..
INCL = -I. -I/opt/sw/gsl-1.10/include -I/opt/sw/vtk-5.0.3/include/vtk-5.0 -I/usr/include/vtk-5.0 -I$(HOME)/include/vtk-5.0 -I$(HOME)/include -I$(CODEDIR)
CPPFLAGS = -Wno-deprecated $(INCL) -O3 -g
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkDICOMParser -lvtkpng -lvtkjpeg -lvtktiff -lvtksys -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
VTKLIBS = -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -ldl
FFTWLIBS = -lfftw3 -lm
GSLLIBS = -lgsl -lgslcblas
LIBS = -L/opt/sw/gsl-1.10/lib -L/opt/sw/vtk-5.0.3/lib $(LDFLAGS) $(VTKLIBS) $(GSLLIBS)  


%.o : $(CODEDIR)/%.cpp
	$(CPPC) $(CPPFLAGS) -c $< -o $@

%.o : %.cpp
	$(CPPC) $(CPPFLAGS) -c $< -o $@

NONSPH = io_data_vtk.o dataLL.o vect.o ioGlobals.o ioDomain.o misc.o
SIM = customSim.o customOutput.o
COMPRESS = $(NONSPH) sphCompress.o $(SIM)
INCOMPRESS = $(NONSPH) sphIncompress.o $(SIM)

run: $(INCOMPRESS) run.o setup createPVD
	$(CPPC) $(CPPFLAGS) -o run $(INCOMPRESS) run.o $(LIBS)
customSim.o: customSim.cpp 
	$(CPPC) $(CPPFLAGS) -c customSim.cpp -o customSim.o 
setup: $(INCOMPRESS) setup.o
	$(CPPC) $(CPPFLAGS) $(LIBS) -o setup $(INCOMPRESS) setup.o $(LIBS)
post: $(INCOMPRESS) post.o
	$(CPPC) $(CPPFLAGS) -o post $(INCOMPRESS) post.o $(LIBS)
convertCSIRO: $(NONSPH) convertCSIRO.o
	$(CPPC) $(CPPFLAGS) -o convertCSIRO $(NONSPH) convertCSIRO.o $(LIBS)
setup.o: setup.cpp 
	$(CPPC) $(CPPFLAGS) -c setup.cpp -o setup.o
createPVD: createPVD.o
	$(CPPC) $(CPPFLAGS) $(LIBS) -o createPVD createPVD.o

clean:
	rm *.o run setup createPVD post convertCSIRO
