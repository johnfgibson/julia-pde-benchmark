function ksbenchmark(Nx)
% ksbenchmark: run a Kuramoto-Sivashinky simulation, benchmark, and plot
%   Nx = number of gridpoints

  Lx = Nx/16*pi;     % spatial domain [0, L] periodic
  dt = 1/16;         % discrete time step 
  T  = 200;          % integrate from t=0 to t=T
  Nt = floor(T/dt);  % total number of timesteps

  x = Lx*(0:Nx-1)/Nx;
  u0 = cos(x) + 0.1*sin(x/8) + 0.01*cos((2*pi/Lx)*x);

  Nruns = 5;
  skip = 1;
  avgtime = 0;

  for r=1:Nruns;
    tic();
    u = ksintegrate(u0, Lx, dt, Nt, nsave);
    cputime = toc()
    if r > skip
      avgtime = avgtime + cputime;
    end
  end

  avgtime = avgtime/(Nruns-skip)
end


function n = ksnorm(u)
% ksnorm: compute the 2-norm of u(x) = sqrt(1/Lx int_0^Lx |u|^2 dx) 
  n = sqrt((u * u') /length(u));
end

function u = ksintegrate(u, Lx, dt, Nt, nsave)
% ksintegrate: integrate kuramoto-sivashinsky equation 
%        u_t = -u*u_x - u_xx - u_xxxx, domain x in [0,Lx], periodic BCs 
%
% inputs
%          u = initial condition (vector of u(x) values on uniform gridpoints))
%         Lx = domain length
%         dt = time step
%         Nt = number of integration timesteps
%      nsave = save every nsave-th time step
%
% outputs
%
%          t = vector of time values (t0, t1, t2, ...) for saved data
%          U = matrix of u(x,t); rows are constant time, cols constant x

  Nx = length(u);                % number of gridpoints
  kx = [0:Nx/2-1 0 -Nx/2+1:-1];  % integer wavenumbers: exp(2*pi*kx*x/L)
  alpha = 2*pi*kx/Lx;            % real wavenumbers:    exp(alpha*x)
  D = i*alpha;                   % D = d/dx operator in Fourier space
  L = alpha.^2 - alpha.^4;       % linear operator -D^2 - D^3 in Fourier space
  G = -0.5*D;                    % -1/2 D operator in Fourier space

  % Express PDE as u_t = Lu + N(u), L is linear part, N nonlinear part.
  % Then Crank-Nicolson Adams-Bashforth discretization is 
  % 
  % (I - dt/2 L) u^{n+1} = (I + dt/2 L) u^n + 3dt/2 N^n - dt/2 N^{n-1}
  %
  % let A = (I - dt/2 L) 
  %     B = (I + dt/2 L), then the CNAB timestep formula is 
  % 
  % u^{n+1} = A^{-1} (B u^n + 3dt/2 N^n - dt/2 N^{n-1})

  % some convenience variables
  dt2  = dt/2;
  dt32 = 3*dt/2;
  A_inv = (ones(1,Nx) - dt2*L).^(-1);
  B     =  ones(1,Nx) + dt2*L;

  Nn  = G.*fft(u.*u); % compute -1/2 d/dx u^2 (spectral), notation Nn = N^n = N(u(n dt))
  Nn1 = Nn;           %                          notation Nn1 = N^{n-1} = N(u((n-1) dt))
  u = fft(u);         % transform u (spectral)
  
  % timestepping loop  
  for n = 1:Nt
      
    Nn1 = Nn;                       % shift N(u) in time: N^{n-1} <- N^n
    Nn  = G.*fft(real(ifft(u)).^2); % compute Nn = N(u) = -1/2 d/dx u^2

    u = A_inv .* (B .* u + dt32*Nn - dt2*Nn1);

  end
  u = real(ifft(u))
end

