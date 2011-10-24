
global re0 = 0.1 
global porosity = 1.0

function f=phi(re)
global re0;
global porosity;
A = 18*re0;
remid = A/18;
beta = 3.7 - 0.65*exp(-((1.5-log10(re)).^2)/2.0);
factor = porosity.^(-beta);
f = 0.3969*re.^2 + 6.048*re.^1.5 + 23.04*re - (4/3)*A./factor;
endfunction

[re, fval, info] = fsolve (@phi, re0)
d = 0.025*2;
visc = 0.0218375/sqrt(re0)
v = re*visc/(d*porosity)
