# julia-pde-benchmark
Benchmarking a simple PDE integration algorithm in Julia and other languages

We benchmark a very simple numerical integration algorithm for a 1d nonlinear partial differential equation 
(PDE) in Julia, Python, Matlab, C++, C, and Fortran. Results include execution time versus size of discretized 
system, and execution time versus lines of code. 

The PDE is the Kuramoto-Sivashinky (KS) equation is the 1d nonlinear PDE  

u_t + u_xx + u_xxxx + u u_x = 0

where x is space and t is time, and subscripts indicate differentiation. The algorithm uses Fourier decomposition in space and 2nd-order Crank-Nicolson, Adams-Bashforth semi-implicit finite-differencing in time, with collocation computation of the nonlinear term. Implementations in all languages utilize the same FFTW library for fast Fourier transforms, so the benchmark is meant to compare the language-specific overheads for things like index bounds checking and allocation of temporary arrays. 

The results show that Julia is competitive with Python and Matlab in line count (about twenty lines of code) 
and C++ and Fortran for execution speed. 

Complete results and discussion are in [1-Kuramoto-Sivashinksy-benchmark.ipynb](1-Kuramoto-Sivashinksy-benchmark.ipynb)

Funding for this research was provided by the National Science Foundation CBET award number 1554149. 
