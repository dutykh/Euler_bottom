%%% -------------------------------------------------- %%%
%%% Initial condition for the full Euler equations     %%%
%%% -------------------------------------------------- %%%
%%% Last modified: 09/02/2016                          %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

% In the current version we set a localized wave-packet
function v = InitCond ()

    global a0 g h1 tol xi L L1 N Nvec

    k0   = 0.05;      % initial wave number
    D0   = 232.0;     % wavepacket length
    p0   = 0.5*L1;    % wavepacket initial position
    om   = sqrt(g*k0*tanh(k0*h1));
    err  = inf;
    x0   = xi; 
    FS   = @(x) a0*cos(k0*(x + L)).*sech((L + x - p0)/D0); % Free surface elevation
    eta0 = FS(x0(Nvec));
    while (err > tol)
        eta  = [eta0; fliplr(eta0)];
        xb   = x_bot(eta, b(x0));
        xf   = x(eta, b(xb));
        eta1 = FS(xf(Nvec));
        err  = norm(eta1 - eta0, inf);
        x0   = xf;
        eta0 = eta1;
    end % while ()
    phi        = (g/om)*a0*sin(k0*(x0(Nvec) + L)).*sech((L + x0(Nvec) - p0)/D0).*cosh(k0*(h1 + eta0))*sech(k0*h1);
    v          = zeros(2*N, 1);
    v(1:N)     = eta;
    v(N+1:end) = [phi; fliplr(phi)];

end % InitCond ()