set xlabel "Temperature"
set ylabel "Average Energy"
set key right bottom

plot '15.dat' title "15 acids" pt 7, \
	'30.dat' title "30 acids" pt 5, \
	'100.dat' title "100 acids" pt 9