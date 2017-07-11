# julia-pde-benchmark
Benchmarking a simple PDE integration algorithm in Julia and other languages

The Kuramoto-Sivashinky (KS) equation is the 1d nonlinear PDE

  u_t + u_xx + u_xxxx + u u_x = 0

where x is space and t is time, and subscripts indicate differentiation. We choose a periodic domain [0, L_x]
and some initial condition u(x,0). We can form a simple numerical intgration scheme for the KS equation using Fourier decomposition in space and 2nd-order Crank-Nicolson, Adams-Bashforth semi-implicit finite-differencing in time, with collocation computation of the nonlinear term ￼. Here we implement the same algorithm in Python, Matlab, C++, Fortran, and Julia. The codes and a detailed description of the algorithm is given below.
