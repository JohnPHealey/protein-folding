f(x) = a*exp(-x)+c
FIT_LIMIT = 1e-6
fit f(x) 'EvsT.dat' via a
plot 'EvsT.dat', f(x)