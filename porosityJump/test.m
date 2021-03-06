A = 18;
porosity = 0.8;
remid = A/18;
range = 1;
re = remid-range:0.01:remid+range;
beta = 3.7 - 0.65*exp(-((1.5-log10(re)).^2)/2.0);
factor = porosity.^(-beta);
f = 0.3969*re.^2 + 6.048*re.^1.5 + 23.04*re - (4/3)*A./factor;
plot(re,f);
axis([remid-range,remid+range,-1.0,1.0]);
pause;
