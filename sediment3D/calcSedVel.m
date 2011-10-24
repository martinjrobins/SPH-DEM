
global re0 = 0.2 
global porosity = 0.8
global d = 100.0e-06;
global visc = 8.9e-07;
dens_dem = 2500.0
dens_fluid = 1000.0;
global A = d^3*(dens_dem/dens_fluid - 1.0)*9.81/(visc^2);

function f=phi(re)
global porosity;
global A;
remid = A/18.0;
beta = 3.7 - 0.65*exp(-((1.5-log10(re)).^2)/2.0);
factor = porosity.^(-beta);
f = 0.3969*re.^2 + 6.048*re.^1.5 + 23.04*re - (4.0/3.0)*A./factor;
endfunction

[re, fval, info] = fsolve (@phi, re0)
v = re*visc/(d)
