set xlabel "Monte Carlo Time"
set ylabel "Energy"
set key right top

plot 'lowT1.dat' title "Simulation 1" with lines, \
	'lowT2.dat' title "Simulation 2" with lines