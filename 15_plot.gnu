f(x) = a*log(b*x)+c
FIT_LIMIT = 1e-6
fit f(x) '15.dat' via a, b, c
plot '15.dat', f(x)