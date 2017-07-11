# julia-pde-benchmark
Benchmarking a simple PDE integration algorithm in Julia and other languages

The Kuramoto-Sivashinky equation is the simple 1d nonlinear PDE

\begin{equation*}
u_t = -u_{xx} -u_{xxxx} - u u_x
\end{equation*}

where $x$ is space and $t$ is time, and subscripts indicate differentiation. We choose a periodic domain $x \in [0, L_x]$
and some initial condition $u(x,0) = u_0(x)$.
