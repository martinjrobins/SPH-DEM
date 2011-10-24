set term postscript colour 20
set output "plotvel.eps"
set xlabel "time (s)"
set ylabel "z velocity (m/s)"
set key bottom left
plot [10.5:14] "demVel.dat" using 1:4 t "dem velocity","demVel.dat" using 1:4 lw 3 smooth unique t "mean dem velocity","demVel.dat" using 1:($4-$23) lw 3 smooth unique t "mean slip velocity", -0.024476 lw 3 t "expected mean"
