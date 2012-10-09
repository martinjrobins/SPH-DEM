CPPC = mpicxx
CODEDIR = ..
INCL = -I$(HOME)/usr/bin    -I. -I/opt/sw/gsl-1.10/include -I/opt/sw/vtk-5.0.3/include/vtk-5.0 -I/usr/include/vtk-5.2 -I/usr/include/vtk-5.8 -I$(HOME)/usr/include/vtk-5.8 -I$(HOME)/usr/include -I$(CODEDIR)
CPPFLAGS =  -Wno-deprecated $(INCL) -O3
#CPPFLAGS =  -Wno-deprecated $(INCL) -g
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkDICOMParser -lvtkpng -lvtkjpeg -lvtktiff -lvtksys -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
VTKLIBS =  -lexpat -lxml2  -lvtkzlib -lvtklibxml2 -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -ldl -lvtksys
FFTWLIBS = -lfftw3 -lm
GSLLIBS = -lgsl -lgslcblas
LIBS = -L$(HOME)/usr/lib -L$(HOME)/usr/lib/vtk-5.8 -L/usr/mpi/gcc/openmpi-1.4.3/lib64/ -L/opt/sw/gsl-1.10/lib -L/opt/sw/vtk-5.0.3/lib $(LDFLAGS) $(VTKLIBS) $(GSLLIBS)


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
calcPost: $(INCOMPRESS) calcPost.o
	$(CPPC) $(CPPFLAGS) -o calcPost $(INCOMPRESS) calcPost.o $(LIBS)
setup.o: setup.cpp 
	$(CPPC) $(CPPFLAGS) -c setup.cpp -o setup.o
createPVD: createPVD.o
	$(CPPC) $(CPPFLAGS) $(LIBS) -o createPVD createPVD.o

clean:
	rm *.o run setup createPVD post calcPost
