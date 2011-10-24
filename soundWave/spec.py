#! /usr/bin/env python
from scipy.ndimage import *
from pylab import *
import matplotlib.axes3d as p3
from scipy import loadtxt,absolute,isinf,integer,resize,real,imag;
from scipy.fftpack import *

pres = 100;
fml = zeros((pres),dtype=float);
fmldpdt = zeros((pres),dtype=float);
fmldpold = zeros((pres),dtype=float);
figure();
subplot(2,1,1);
hold(False);
subplot(2,1,2);
hold(False);
for i in range(1,322):
   print i

   vin = loadtxt('resultsVGrid0000%03d.dat'%i);
   v = zeros((pres,pres),dtype=float);
   v = vin[:pres,:pres];
   #v = 0.5*v*v*1000*((2.0/300.0)**2);

   ml = v[pres/2,:];
   fmlfft = fft(ml);
   fmlfftr = real(fmlfft);
   fmlffti = imag(fmlfft);

   fml = sqrt(fmlfftr*fmlfftr + fmlffti*fmlffti);
   fmlp = arctan2(fmlffti,fmlfftr);
   fmldpdt = fmlp-fmldpold;
   fmldpdt[fmldpdt>pi/2.0] = fmldpdt[fmldpdt>pi/2.0]-pi;
   fmldold = fmlp;

   fl = arange(0,fml.shape[0],dtype=float);
   subplot(2,1,1);
   plot(fl,fml);
   ax = list(axis());
   ax[1] = pres/2;
   axis(ax);
   subplot(2,1,2);
   plot(fl,fmldpdt);
   ax = list(axis());
   ax[1] = 10;
   axis(ax);
   savefig('resultsSpec0000%03d.png'%i);


