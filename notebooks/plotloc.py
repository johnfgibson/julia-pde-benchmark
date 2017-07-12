from pylab import *

D = loadtxt('timeloc.asc',comments='#');

plot(D[0,1], D[0,0], 'mo', label='Python')
plot(D[1,1], D[1,0], 'yo', label='Julia, naive')
plot(D[2,1], D[2,0], 'bo', label='Matlab')
plot(D[3,1], D[3,0], 'ro', label='Julia, in place')
plot(D[4,1], D[4,0], 'co', label='C++')

xlim(0, 100)
ylim(0, 50)
xlabel('lines of code');
ylabel('cpu time');
#legend(loc='upper left')
legend(loc='upper right')
title('Kuramoto-Shivashinksy simulation')
