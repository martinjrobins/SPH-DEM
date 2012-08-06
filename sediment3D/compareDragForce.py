from numpy import *
from pylab import *
#epsilon = arange(0.37,1,0.001)
epsilon = arange(0.8,0.95,0.05)
dens = 1000.0
visc = 8.9e-4/dens
Rep = 0.85
#dens = 1.1839
#visc = 18.6e-6/dens
#Rep = 3.19 
#dens = 1150.0
#visc = 8.9e-3/dens
Rep = 0.011
dem_dens = 2.5*dens
dem_d = 100.0e-6
dv = Rep*visc/dem_d
if (Rep>1000):
   print "not valid!!!"
Cd = 24.0*(1.0+0.15*Rep**0.687)/Rep
gamma = 3.7-0.65*exp(-(1.5-log10(Rep))**2/2.0)
voidage = epsilon**(-gamma+2)
fStokes = 3*pi*visc*dens*dem_d*dv*epsilon/epsilon
fDiFelice = voidage*Cd*pi*dens*dem_d**2*dv**2/8.0
beta = copy(epsilon)
beta[epsilon<=0.8] = 150.0*((1.0-epsilon[epsilon<=0.8])**2/epsilon[epsilon<=0.8])*(visc*dens/dem_d**2) + 1.75*(dens/dem_d)*dv*(1.0-epsilon[epsilon<=0.8])/epsilon[epsilon<=0.8] 
beta[epsilon>0.8] = (3.0/4.0)*(Cd/dem_d**2)*dens*visc*(1.0-epsilon[epsilon>0.8])*epsilon[epsilon>0.8]**(-2.7)*Rep
fErgun = (4.0/3.0)*pi*(dem_d/2.0)**3*beta*dv/(1.0-epsilon)
krRep = dv*dem_d/visc
print krRep
fKhanRichardson = pi*(dem_d/2.0)**2*dens*dv**2*(1.84*krRep**(-0.31)+0.293*krRep**0.06)**3.45*epsilon/epsilon

print fDiFelice
plot(epsilon,fStokes,label='Stokes')
plot(epsilon,fDiFelice,label='DiFelice')
plot(epsilon,fKhanRichardson,label='Khan and Richardson')
plot(epsilon,fErgun,label='Ergun + Wen and Yu')
xlabel("porosity")
ylabel("drag force")
yscale('log')
legend(loc=4)



epsilon = 1.0 
dem_dens = 2.5*dens
dem_d = 100.0e-6
Rep = arange(0.01,10,0.01)
visc = 8.9e-4/dens
dv = Rep*visc/dem_d
Cd = 24.0*(1.0+0.15*Rep**0.687)/Rep
gamma = 3.7-0.65*exp(-(1.5-log10(Rep))**2/2.0)
voidage = epsilon**(-gamma+2)
fStokes = 3*pi*visc*dens*dem_d*dv
fDiFelice = voidage*Cd*pi*dens*dem_d**2*dv**2/8.0
beta = (3.0/4.0)*(Cd/dem_d**2)*dens*visc*(1.0-epsilon)*epsilon**(-2.7)*Rep
fErgun = (4.0/3.0)*pi*(dem_d/2.0)**3*beta*dv
krRep = Rep
fKhanRichardson = pi*(dem_d/2.0)**2*dens*dv**2*(1.84*krRep**(-0.31)+0.293*krRep**0.06)**3.45
figure(2)
plot(Rep,fStokes,label='Stokes')
plot(Rep,fDiFelice,label='DiFelice')
plot(Rep,fKhanRichardson,label='Khan and Richardson')
plot(Rep,fErgun,label='Ergun + Wen and Yu')
xlabel("Reynolds Number")
ylabel("drag force")
yscale('log')
legend(loc=4)

figure(3)
plot(Rep,fStokes,label='Stokes')
plot(Rep,fDiFelice,label='DiFelice')
plot(Rep,fKhanRichardson,label='Khan and Richardson')
plot(Rep,fErgun,label='Ergun + Wen and Yu')
xlabel("Reynolds Number")
ylabel("drag force")
legend(loc=2)
show()
