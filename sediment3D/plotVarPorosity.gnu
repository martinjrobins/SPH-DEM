set term postscript colour 
set size 0.6,0.7
set xlabel "time (s)"
set ylabel "z velocity (m/s)"
set output "varPorosity0.01.eps"
plot [9.0:16] "multi0.7gridBlock10_0.01/demVel.dat" using 1:($4+0.7*0.017132)/(0.7*0.017132) smooth unique t "0.7", "multi0.8gridBlock10_0.01/demVel.dat" using 1:($4+0.8*0.024476)/(0.8*0.024476) smooth unique t "0.8", "multi0.9gridBlock10_0.01/demVel.dat" using 1:($4+0.9*0.033499)/(0.9*0.033499) smooth unique t "0.9", "multi0.95gridBlock10_0.01/demVel.dat" using 1:($4+0.95*0.038676)/(0.95*0.038676) smooth unique t "0.95"
set output "varPorosity1.0.eps"
plot [9.0:16] "multi0.7gridBlock10_1.0/demVel.dat" using 1:($4+0.7*0.15685)/(0.7*0.15685) smooth unique t "0.7", "multi0.8gridBlock10_1.0/demVel.dat" using 1:($4+0.21754*0.8)/(0.8*0.21754) smooth unique t "0.8", "multi0.9gridBlock10_1.0/demVel.dat" using 1:($4+0.9*0.28681)/(0.9*0.28681) smooth unique t "0.9", "multi0.95gridBlock10_1.0/demVel.dat" using 1:($4+0.95*0.32412)/(0.95*0.32412) smooth unique t "0.95"
set output "varPorosityTest.eps"
plot [9.0:16] "multi0.7gridBlock10_0.01/demVel.dat" using 1:4 smooth unique t "0.7", "multi0.8gridBlock10_0.01/demVel.dat" using 1:4 smooth unique t "0.8", "multi0.9gridBlock10_0.01/demVel.dat" using 1:4 smooth unique t "0.9", "multi0.95gridBlock10_0.01/demVel.dat" using 1:4 smooth unique t "0.95"
