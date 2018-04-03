set xlabel "Temperature"
set ylabel "Average Energy"

f(x) = a*atan(b*(x-2.4))-130
FIT_LIMIT = 1e-8
fit f(x) '100.dat' via a, b
plot '100.dat' title "" pt 7, f(x) title ""