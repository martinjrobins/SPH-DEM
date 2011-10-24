

Csph gsph;
Csimulation gsimulation;
Cdata gdata;
Cio_data gio_data;

void main() {
	//read in particles from file
	the_data = gio_data.read(filename);
	//give them to data structure
	gdata.construct(the_data);
	//start simulation
	gsimulation.run();
	//clean up and get out of here
}
