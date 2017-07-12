from pylab import *

T = loadtxt('timings.asc',comments='%');
N = T[:,0]; # gridsize



#loglog(N, T[:,3], 'g-o', label='Octave')
loglog(N, T[:,4], 'm-o', label='Python')
#loglog(N, T[:,5], 'r-o', label='Chflow')
loglog(N, T[:,2], 'b-o', label='Matlab')
loglog(N, T[:,1], 'c-o', label='C++')
loglog(N, T[:,6], 'y-o', label='Julia, naive')
loglog(N, T[:,8], 'r-o', label='Julia, in place')
#loglog(N, T[:,7], 'k-o', label='Julia, unrolled')

loglog(N, 0.0000125*N*log10(N), 'k--', label='N log N')
ylim(1e-03, 1e02)
xlim(10, 1e06)
xlabel('N');
ylabel('cpu time');
#legend(loc='upper left')
legend(loc='lower right')
title('Kuramoto-Shivashinksy simulation, N gridpoints')
