
re = 0.1 
porosity = 1.0
dens = 1000
densp = 2500

visc = 8.92635148e-4;
beta = 3.7 - 0.65*exp(-((1.5-log10(re)).^2)/2.0);
factor = porosity.^(-beta);
A =  (3/4)*factor*(0.3969*re.^2 + 6.048*re.^1.5 + 23.04*re);

d = (A*(visc^2)/(dens*(densp-dens)*9.81))^(1.0/3.0)
