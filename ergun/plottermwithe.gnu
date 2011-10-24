set term postscript
set output "test2.eps"
plot "-" using 1:($2*$1)  w lp ps 3 lw 3 t "real velocity times porosity"
#  X     Y     
0.6 0.011337
0.7 0.017132
0.8 0.024476
0.9 0.033499
0.99 0.043150
end


