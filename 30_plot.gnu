set xlabel "Temperature"
set ylabel "Average Energy"

f(x) = a*atan(b*(x-1.2))-70
FIT_LIMIT = 1e-8
fit f(x) '30.dat' via a, b
plot '30.dat' title "" pt 7, f(x) title ""