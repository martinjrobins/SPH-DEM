CPPC = g++
CPPFLAGS = -Wno-deprecated -O3 -I. -I/usr/include/vtk-5.0
HDF5DIR = /usr/local/hdf5
H5PARTDIR = /usr/local/h5part
VTKLIBS = -lvtkCommon -lvtkIO


%.o : %.cpp
	$(CPPC) $(CPPFLAGS) -c $< -o $@

TESTDATAOBJ = vect.o particle.o testData.o

NONSPH = data.o vect.o io_data_vtk.o setupSod.o
T1DSim = data.o vect.o sphCompress.o io_data_vtk.o 1DSim.o
COMPRESS = $(NONSPH) sphCompress.o run.o
INCOMPRESS = $(NONSPH) sphIncompress.o run.o

couette/run: $(INCOMPRESS) couette/setupCouette
	$(CPPC) $(CPPFLAGS) -o couette/run $(INCOMPRESS) $(VTKLIBS) -I. -Icouette
couette/setupCouette: $(NONSPH)
	$(CPPC) $(CPPFLAGS) -o couette/setupCouette $(NONSPH) setupCouette.o $(VTKLIBS) -I. -Icouette
T1DSim: $(T1DSim)
	$(CPPC) $(CPPFLAGS) -o 1DSim $(T1DSim) -lvtkCommon -lvtkIO
setupSod: $(SETUPSOD)
	$(CPPC) $(CPPFLAGS) -o setupSod $(SETUPSOD) -lvtkCommon -lvtkIO
testData: $(TESTDATAOBJ) 
	g++ -o testData $(TESTDATAOBJ) 
testData.o: testData.cpp 
	$(CPPC) $(CPPFLAGS) -c testData.cpp -o testData.o
io_data_vtk.o: io_data_vtk.cpp
	$(CPPC) $(CPPFLAGS) -c io_data_vtk.cpp -o io_data_vtk.o 
io_data_h5part.o: io_data_h5part.cpp
	$(CPPC) $(CPPFLAGS) -c io_data_h5part.cpp -o io_data_h5part.o -I $(HDF5DIR)/include/ -I $(H5PARTDIR)/include/
clean:
	rm *.o
