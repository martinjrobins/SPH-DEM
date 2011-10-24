from numpy import *
from pylab import *
epsilon = arange(0.5,1,0.001)
dv = 1.0
dens = 1000
visc = 8.9e-4/dens
dem_dens = 2.5*dens
dem_d = 100e-6
Rep = dem_d*dv*epsilon/visc
if any(Rep>1000):
   print "not valid!!!"
Cd = 24.0*(1.0+0.15*Rep**0.687)/Rep
gamma = 3.7-0.65*exp(-(1.5-log10(Rep))**2/2.0)
voidage = epsilon**(-gamma+2)
fStokes = 3*pi*visc*dens*dem_d*dv*epsilon/epsilon
fDiFelice = voidage*Cd*pi*dens*dem_d**2*dv**2/8.0
beta = copy(epsilon)
beta[epsilon<=0.8] = 150*((1-epsilon)/epsilon)*(visc*dens/dem_d**2) + 1.75*(dens/dem_d)*dv 
beta[epsilon>0.8] = (3.0/4.0)*(Cd/dem_d**2)*dens*visc*epsilon**(-2.7)*Rep
fErgun = (4.0/3.0)*pi*(dem_d/2.0)**3*beta*dv
krRep = dv*dem_d/visc
fKhanRichardson = pi*(dem_d/2.0)**2*dens*dv**2*(1.84*krRep**(-0.31)+0.293*krRep**0.06)**3.45*epsilon/epsilon
plot(epsilon,fStokes,label='Stokes')
plot(epsilon,fDiFelice,label='DiFelice')
plot(epsilon,fKhanRichardson,label='Khan and Richardson')
plot(epsilon,fErgun,label='Ergun + Wen and Yu')
xlabel("porosity")
ylabel("drag force")
legend()


epsilon = 1.0 
dv = 1.0
dens = 1000
dem_dens = 2.5*dens
dem_d = 100e-6
Rep = arange(0.01,10,0.01)
visc = dv*dem_d*epsilon/Rep
Cd = 24.0*(1.0+0.15*Rep**0.687)/Rep
CdTL = 24.0/Rep + 4.152/Rep**0.343 + 0.413/(1+16300*Rep**(-1.09))
gamma = 3.7-0.65*exp(-(1.5-log10(Rep))**2/2.0)
voidage = epsilon**(-gamma+2)
fStokes = 3*pi*visc*dens*dem_d*dv
fDiFelice = voidage*Cd*pi*dens*dem_d**2*dv**2/8.0
fDiFeliceTL = voidage*CdTL*pi*dens*dem_d**2*dv**2/8.0
beta = (3.0/4.0)*(Cd/dem_d**2)*dens*visc*epsilon**(-2.7)*Rep
fErgun = (4.0/3.0)*pi*(dem_d/2.0)**3*beta*dv
betaTL = (3.0/4.0)*(CdTL/dem_d**2)*dens*visc*epsilon**(-2.7)*Rep
fErgunTL = (4.0/3.0)*pi*(dem_d/2.0)**3*betaTL*dv
krRep = Rep
fKhanRichardson = pi*(dem_d/2.0)**2*dens*dv**2*(1.84*krRep**(-0.31)+0.293*krRep**0.06)**3.45
figure(2)
plot(Rep,fStokes,label='Stokes')
plot(Rep,fDiFelice,label='DiFelice')
plot(Rep,fDiFeliceTL,label='DiFelice TL')
plot(Rep,fKhanRichardson,label='Khan and Richardson')
plot(Rep,fErgun,label='Ergun + Wen and Yu')
plot(Rep,fErgunTL,label='Ergun + Wen and Yu TL')
xlabel("Reynolds Number")
ylabel("drag force")
yscale('log')
legend()

figure(3)
plot(Rep,fDiFelice,label='DiFelice')
plot(Rep,fDiFeliceTL,label='DiFelice TL')
plot(Rep,fKhanRichardson,label='Khan and Richardson')
plot(Rep,fErgun,label='Ergun + Wen and Yu')
plot(Rep,fErgunTL,label='Ergun + Wen and Yu TL')
xlabel("Reynolds Number")
ylabel("drag force")
yscale('log')
xscale('log')
legend()
show()
