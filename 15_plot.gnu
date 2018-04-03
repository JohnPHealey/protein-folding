set xlabel "Temperature"
set ylabel "Average Energy"

f(x) = a*atan(b*(x-2))-25
FIT_LIMIT = 1e-6
fit f(x) '15.dat' via a, b
plot '15.dat' title "" pt 7, f(x) title ""