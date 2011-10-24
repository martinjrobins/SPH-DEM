set term postscript colour 
set size 0.6,0.7
set xlabel "time (s)"
set ylabel "z velocity (m/s)"
p=0.8
es = 0.024476
set output "plotBlock0.01.eps"
plot [9.0:16] "multi0.8gridBlock10_0.01/demVel.dat" using 1:4 smooth unique t "10", "multi0.8gridBlock15_0.01/demVel.dat" using 1:4 smooth unique t "15", "multi0.8gridBlock20_0.01/demVel.dat" using 1:4 smooth unique t "20", "multi0.8gridBlock25_0.01/demVel.dat" using 1:4 smooth unique t "25", -p*0.024476 lw 3 t "expected slip"
set output "plotBlock1.00.eps"
plot [9.0:16] "multi0.8gridBlock10_1.0/demVel.dat" using 1:4 smooth unique t "10", "multi0.8gridBlock15_1.0/demVel.dat" using 1:4 smooth unique t "15", "multi0.8gridBlock20_1.0/demVel.dat" using 1:4 smooth unique t "20", "multi0.8gridBlock25_1.0/demVel.dat" using 1:4 smooth unique t "25", -p*0.21754 lw 3 t "expected slip"
set output "plotBlockError0.01.eps"
plot [11.0:16] "multi0.8gridBlock10_0.01/demVel.dat" using 1:($4+p*es)/(p*es) smooth unique t "10", "multi0.8gridBlock15_0.01/demVel.dat" using 1:($4+p*es)/(p*es) smooth unique t "15", "multi0.8gridBlock20_0.01/demVel.dat" using 1:($4+p*es)/(p*es) smooth unique t "20", "multi0.8gridBlock25_0.01/demVel.dat" using 1:($4+p*es)/(p*es) smooth unique t "25"
set output "plotBlockError1.00.eps"
plot [11:16] "multi0.8gridBlock10_1.0/demVel.dat" using 1:($4+p*es)/(p*es) smooth unique t "10", "multi0.8gridBlock15_1.0/demVel.dat" using 1:($4+p*es)/(p*es) smooth unique t "15", "multi0.8gridBlock20_1.0/demVel.dat" using 1:($4+p*es)/(p*es) smooth unique t "20", "multi0.8gridBlock25_1.0/demVel.dat" using 1:($4+p*es)/(p*es) smooth unique t "25"

set xlabel "h / d"
set ylabel "% term. vel. error"
set output "error.eps"
plot 'error.dat' using (1/$1*1.5/0.05):($2*100) w lp notitle
