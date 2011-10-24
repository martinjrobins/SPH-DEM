b = 6.0*pi*0.0218375*1000.0*0.025
v = (4.0/3.0)*pi*0.025**3
m = v*8000
set term postscript colour 20
set xlabel "time (s)"
set ylabel "Vertical Velocity"
set parametric
set trange [10:13]

set output "termVelError.eps"
plot [10:13][11:11.4] "single8/demVel.dat" using 1:(($4+(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m)))/(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m))) w p ps 1.5 lw 2 title "8","single10/demVel.dat" using 1:(($4+(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m)))/(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m))) w p ps 1.5 lw 2 title "10","single15/demVel.dat" using 1:(($4+(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m)))/(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m))) w p ps 1.5 lw 2 title "15","single20/demVel.dat" using 1:(($4+(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m)))/(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m))) w p ps 1.5 lw 2 title "20","single25/demVel.dat" using 1:(($4+(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m)))/(7000.0*v*9.81/b)*(1.0-exp(-b*($1-10)/m))) w p ps 1.5 lw 2 title "25"


set output "termVelErrorPlot.eps"
set xrange [9:26]
set xlabel "Resolution"
set ylabel "% Error"
plot "-" using 1:(100*$2)  w lp ps 3 lw 3 notitle
#  X     Y     
10 0.04
15 0.1
20 0.15
25 0.18
end

set output "termVelTwoWay.eps"
set xlabel "time (s)"
set ylabel "Vertical Velocity"
set key bottom left
set size 1,1
set origin 0,0
set multiplot
set size 1,1
set origin 0,0
plot [10:13][9.5:11] "single10/demVel.dat" using 1:4 w p ps 1.5 lw 2 title "10","single15/demVel.dat" using 1:4 w p ps 1.5 lw 2 title "15","single20/demVel.dat" using 1:4 w p ps 1.5 lw 2 title "20","single25/demVel.dat" using 1:4 w p ps 1.5 lw 2 title "25",t,-(7000.0*v*9.81/b)*(1.0-exp(-b*(t-10)/m)) w l lw 2 title "Stokes Law"

set size 0.5,0.5
set origin 0.45,0.45
unset xlabel
set ylabel "Slip Velocity"
set xtics 0.5 
set ytics 0.2 
plot [10:13][9.5:11] "single10/demVel.dat" using 1:($4-$23) w p ps 1 lw 1 notitle,"single15/demVel.dat" using 1:($4-$23) w p ps 1 lw 1 notitle,"single20/demVel.dat" using 1:($4-$23) w p ps 1 lw 1 notitle,"single25/demVel.dat" using 1:($4-$23) w p ps 1lw 1 notitle,t,-(7000.0*v*9.81/b)*(1.0-exp(-b*(t-10)/m)) w l lw 1 notitle 


