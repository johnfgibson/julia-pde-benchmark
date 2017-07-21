```
ksintegrate: integrate kuramoto-sivashinsky equation (Julia)
       u_t = -u*u_x - u_xx - u_xxxx, domain x in [0,Lx], periodic BCs 

 inputs
          u = initial condition (vector of u(x) values on uniform gridpoints))
         Lx = domain length
         dt = time step
         Nt = number of integration timesteps
      nsave = save every nsave-th time step

 outputs

          u = final state, vector of u(x, Nt*dt) at uniform x gridpoints

This implementation has two improvements over ksintegrateNaive.jl. It uses 
   (1) in-place FFTs
   (2) loop fusion: Julia can translate arithmetic vector expressions in dot syntax 
       to single for loop over the components, which should be much faster than
       constructing a temporary vector for each operation in the vector expression. 
   
```

function ksintegrateInplace(u, Lx, dt, Nt);
    u = (1+0im)*u                       # force u to be complex
    Nx = length(u)                      # number of gridpoints
    kx = vcat(0:Nx/2-1, 0:0, -Nx/2+1:-1)# integer wave #s,  exp(2*pi*i*kx*x/L)
    alpha = 2*pi*kx/Lx                  # real wavenumbers, exp(i*alpha*x)
    D = 1im*alpha                       # spectral D = d/dx operator 
    L = alpha.^2 - alpha.^4             # spectral L = -D^2 - D^4 operator
    G = -0.5*D                          # spectral -1/2 D operator, to eval -u u_x = 1/2 d/dx u^2

    # convenience variables
    dt2  = dt/2
    dt32 = 3*dt/2
    A_inv = (ones(Nx) - dt2*L).^(-1)
    B     =  ones(Nx) + dt2*L

    # compute FFTW plans
    FFT! = plan_fft!(u, flags=FFTW.ESTIMATE)
    IFFT! = plan_ifft!(u, flags=FFTW.ESTIMATE)

    # compute nonlinear term Nu == -u u_x and Nuprev (Nu at prev timestep)
    Nu     = G.*fft(u.^2); # Nu == -1/2 d/dx (u^2) = -u u_x
    Nuprev = copy(Nu);      # use Nuprev = Nu at first time step
    FFT!*u;

    # timestepping loop
    for n = 0:Nt

        Nuprev .= Nu   # shift nonlinear term in time
        Nu .= u         # put u into N in prep for comp of nonlineat
        
        IFFT!*Nu;       # transform Nu to gridpt values, in place
        Nu .= Nu.*Nu;   # collocation calculation of u^2
        FFT!*Nu;        # transform Nu back to spectral coeffs, in place

        Nu .= G.*Nu;

        # loop fusion! Julia translates the folling line of code to a single for loop. 
        u .= A_inv .*(B .* u .+ dt32.*Nu .- dt2.*Nuprev); 
    end

    IFFT!*u
    real(u)
end
