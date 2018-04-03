set xlabel "Monte Carlo Time"
set ylabel "Energy"
set key right top

plot 'Esame.dat' title "E_{i,j} = 1" with lines, \
	'Ediff.dat' title "E_{i,j} = +/- 1" with lines