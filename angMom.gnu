set xlabel "Time (s)"
set ylabel "Normalised Angular Momentum"
plot 'resultsGlobals0.dat' using 3:($6/sqrt($7*4000*3/4)) with l notitle
