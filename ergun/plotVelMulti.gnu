set term postscript colour 20
set output "plotvel.eps"
set xlabel "time (s)"
set ylabel "z velocity (m/s)"
p=0.8
set key top left
plot [10.5:16] "demVel.dat" using 1:4 t "dem velocity","demVel.dat" using 1:4 lw 3 smooth unique t "mean dem velocity","demVel.dat" using 1:($4-$23)*p lw 3 smooth unique t "mean slip velocity", -p*0.024476 lw 3 t "expected slip"


set key bottom left
set output "plotvelzoom.eps"
plot [10.5:14][-0.05:0.0] "demVel.dat" using 1:4 t "dem velocity","demVel.dat" using 1:4 lw 3 smooth unique t "mean dem velocity","demVel.dat" using 1:($4-$23)*p lw 3 smooth unique t "mean slip velocity", -p*0.024476 lw 3 t "expected slip"

