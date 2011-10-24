CPPC = mpic++
CODEDIR = ..
INCL = -I. -I/usr/include/vtk-5.0 -I$(HOME)/include/vtk-5.0 -I$(HOME)/include -I$(CODEDIR)
CPPFLAGS = -Wno-deprecated $(INCL) -O3
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
VTKLIBS = -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -ldl
FFTWLIBS = -lfftw3 -lm
GSLLIBS = -lgsl -lgslcblas
LIBS = -L$(HOME)/lib $(VTKLIBS) $(GSLLIBS) $(FFTWLIBS) 


%.o : $(CODEDIR)/%.cpp
	$(CPPC) $(CPPFLAGS) -c $< -o $@

%.o : %.cpp
	$(CPPC) $(CPPFLAGS) -c $< -o $@

NONSPH = io_data_vtk.o dataLL.o vect.o ioGlobals.o ioDomain.o
SIM = customSim.o customOutput.o
COMPRESS = $(NONSPH) sphCompress.o $(SIM)
INCOMPRESS = $(NONSPH) sphIncompress.o $(SIM)

run: $(INCOMPRESS) run.o setup
	$(CPPC) $(CPPFLAGS) -o run $(INCOMPRESS) run.o $(LIBS)
customSim.o: customSim.cpp 
	$(CPPC) $(CPPFLAGS) -c customSim.cpp -o customSim.o 
setup: $(INCOMPRESS) setup.o
	$(CPPC) $(CPPFLAGS) $(LIBS) -o setup $(INCOMPRESS) setup.o $(LIBS)
post: $(NONSPH) customOutput.o post.o
	$(CPPC) $(CPPFLAGS) -o post $(NONSPH) customOutput.o post.o $(LIBS)
setup.o: setup.cpp 
	$(CPPC) $(CPPFLAGS) -c setup.cpp -o setup.o 

clean:
	rm *.o run setup post
