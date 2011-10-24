set terminal postscript monochrome
set output 'taylorGreenMaxV.ps'
set xlabel 'Time (s)'
set ylabel 'Maximum Velocity'
set logscale y
plot 'resultsGlobals0.dat' using 3:20 title 'sph 50x50', exp(-8*3.14*3.14*(x-2.407)/100) title 'Re=100', exp(-8*3.14*3.14*(x-2.407)/93) title 'Re=93'

set terminal x11
set output 
set xlabel 'Time (s)'
set ylabel 'Maximum Velocity'
set logscale y
plot 'resultsGlobals0.dat' using 3:20 title 'sph 50x50', exp(-8*3.14*3.14*(x-2.407)/100) title 'Re=100', exp(-8*3.14*3.14*(x-2.407)/93) title 'Re=93'

