const double H1 = 0.41; //the water height from 0<x<L1 is H1+H2+H3
const double H2 = 0.9;  //the water height from L1<=x<L1+L2 is H1+H2
const double H3 = 0.0;  //etc..
const double L1 = 0.0;  //
const double L2 = 1.0;  //from L1+L2+L3<x<L1+L2+L3+L4 is the start of the fully 3D bit
const double L3 = 20.0;  //from L1+L2+L3+L4<x<L1+L2+L3+L4+L5 is the contraction
const double L4 = 1.0;
const double L5 = 1.0;
const double WIDTH = 1.0;

const int NCPU = 7; //number of cpus to use. Its MPI so go nuts! :)

const double MAXTIME = 7.0;
const double WALLUP = 0.2; //time the gate starts to be removed
const double WALLSPEED = 2.0; 
const int OUTSTEP = 1000;

//#define FLUID_GATE
const int gateWidth = 4; //gate width measured in PSEP (particle separation)

//PY and RY define resolution. the periodic section is 1/RY the width of the
//channel and there are 2*PY particles across this periodic section
const int PY = 4;
const int RY = 0;


#define CREATE_TIP
