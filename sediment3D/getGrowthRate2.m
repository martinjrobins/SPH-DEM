function r = phi(viscosity)
   global porosity;
   global fluid_density;
   global dem_density;
   global L;
   global n;
   k = sqrt(2.0*4.0*pi^2/(L^2));
   rho1 = fluid_density;
   rho2 = fluid_density*porosity+dem_density*(1.0-porosity);
   mu1 = viscosity*rho1;
   mu2 = mu1*(1.0-(1.0-porosity)/0.63)^(-2.5*0.63);
   #mu2 = mu1*(1.0-2.5*(1.0-porosity));
   alpha1 = rho1/(rho1+rho2);
   alpha2 = rho2/(rho1+rho2);
   v1 = mu1/rho1;
   v2 = mu2/rho2;
   q1 = sqrt(k^2+n/v1);
   q2 = sqrt(k^2+n/v2);
   T = 0;
   g = 9.81;
   
   r = -((g*k/n^2)*((alpha1-alpha2)+k^2*T/(g*(rho1+rho2)))+1.0)*(alpha2*q1 + alpha1*q2 - k) - 4.0*k*alpha1*alpha2 + (4.0*k^2/n)*(alpha1*v1 - alpha2*v2)*(alpha2*q1-alpha1*q2+k*(alpha1-alpha2)) + (4.0*k^3/n^2)*(alpha1*v1-alpha2*v2)^2*(q1-k)*(q2-k);
   #r = (k/n^2)*(rho2-rho1)*g-T*k^3/n^2-(2*k^2/n)*(mu2+mu1)-rho2-rho1;
endfunction

global n = 5.08114;
global porosity = 0.93;
viscosity0 = 7.73913043478e-06;
#global viscosity = 8.9e-07;
#global viscosity = 0.02;
#global fluid_density = 1000.0;
global fluid_density = 1150.0;
global dem_density = 2500.0;
#global dem_density = 8000.0;
#global L = 1.0;
global L = 0.004;

#[n,obj,info,iter,nf,lambda] = sqp(n0,@phi)
[viscosity, fval, info] = fsolve (@phi, viscosity0)
viscosity

k = sqrt(2.0*4.0*pi^2/(L^2));
rho1 = fluid_density
rho2 = fluid_density*porosity+dem_density*(1.0-porosity)
mu1 = viscosity*rho1
mu2 = mu1*(1.0-(1.0-porosity)/0.63)^(-2.5*0.63)
mu22 = mu1*((porosity-0.37)/(1-0.37))^(-2.5*(1-0.37))
alpha1 = rho1/(rho1+rho2);
alpha2 = rho2/(rho1+rho2);
v1 = mu1/rho1;
v2 = mu2/rho2;
q1 = sqrt(k^2+n/v1);
q2 = sqrt(k^2+n/v2);
T = 0.0;
g = 9.81;
   

r = -((g*k/n^2)*((alpha1-alpha2)+k^2*T/(g*(rho1+rho2)))+1)*(alpha2*q1 + alpha1*q2 - k) - 4*k*alpha1*alpha2 + (4*k^2/n)*(alpha1*v1 - alpha2*v2)*(alpha2*q1-alpha1*q2+k*(alpha1-alpha2)) + (4*k^3/n^2)*(alpha1*v1-alpha2*v2)^2*(q1-k)*(q2-k)
#r = (k/n^2)*(rho2-rho1)*g-T*k^3/n^2-(2*k^2/n)*(mu2+mu1)-rho2-rho1

A = (rho2-rho1)/(rho1+rho2)
growth = sqrt(A*9.81*k)
