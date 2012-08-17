from scipy import integrate
from scipy import linspace

from scipy.optimize import fsolve
from numpy import exp,log10
from pylab import *
from math import pi

def f(re,ar,porosity):
   beta = 3.7 - 0.65*exp(-((1.5-log10(re))**2.0)/2.0);
   factor = porosity**(1+beta);
   f = 0.3969*re**2 + 6.048*re**1.5 + 23.04*re - (4.0/3.0)*ar*factor;
   return f


def deriv(u,t,porosity,fluid,d,dem_dens):
   x = u[0]; v = u[1]
   re = d*abs(v)/fluid.visc
   if (v==0):
      cd = 0
   else:
      cd = (0.63 + 4.8/sqrt(re))**2*v/v
   beta = 3.7 - 0.65*exp(-((1.5-log10(re))**2.0)/2.0);
   vol = (1.0/6.0)*pi*d**3
   dem_mass = vol*dem_dens
   reduced_mass = vol*porosity*(dem_dens - fluid.dens)
   a = (1.0/4.0)*pi*d**2
   f = -reduced_mass*9.81 + 0.5*fluid.dens*(v**2)*cd*a*(porosity**(-beta))
   print x,v,f
   return (v,f/dem_mass)



class fluid:
   pass


water_glycerol = fluid();
water_glycerol.dens = 1150.0;
water_glycerol.visc = 8.9e-3/water_glycerol.dens;
water_glycerol.name = "water_glycerol";
water = fluid();
water.dens = 1000.0;
water.visc = 8.9e-4/water.dens;
water.name = "water";
water2 = fluid();
water2.dens = 1000.0;
water2.visc = 0.6*8.9e-4/water.dens;
water2.name = "water2Re2Re";
water3 = fluid();
water3.dens = 1000.0;
water3.visc = 0.24*8.9e-4/water.dens;
water3.name = "water10Re";
air = fluid();
air.dens = 1.1839;
air.visc = 18.6e-6/air.dens;
air.name = "air";

d = 100.0e-6;
dem_dens = 2500.0;

for fluid in [water]:
   for porosity in [0.85]:
      ar = d**3*(dem_dens/fluid.dens - 1.0)*9.81/(fluid.visc**2);
      re0 = ar/18.0;
      rePlus = fsolve(f,re0,(ar,porosity));
      re = rePlus;

      t = linspace(0,0.5,num=1000)
      u = integrate.odeint(deriv, (0,0), t, args=(porosity,fluid,d,dem_dens))
      v_2 = t*re*fluid.visc/d/t

      fil = open('analyticalV.dat','w')
      for i in range(1000):
         fil.write(str(t[i])+' '+str(u[i,1])+'\n')
      fil.close()
