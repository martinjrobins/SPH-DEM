f(x) = c-a*exp(x*b)-0.024476*x
fit[10.8:12.2] f(x) 'minDem.dat' using 1:2 via a,b,c

plot [10.8:12.2] 'minDem.dat' using 1:2,f(x) 
