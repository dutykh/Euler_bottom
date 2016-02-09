%%% -------------------------------------------------- %%%
%%% RHS of full Euler equations in a conformal domain  %%%
%%% -------------------------------------------------- %%%
%%% Last modified: 09/02/2016                          %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function rhs = RHS (~, v)

    % Declaration of needed globals:
    global g k kf tol xi N

    rhs     = zeros(2*N,1); % declare the right-hand side
                            % and allocate the memory
    gamma   = v(1:N);       % retrieve the free surface elevation
    phi     = v(N+1:end);   % retrive the velocity potential

    gam_hat = fft(gamma);
    phi_hat = fft(phi);

    % compute iteratively the bathymetry function
    err   = inf;
    h     = b(xi);
    while (err > tol)
        xb  = x_bot(gamma, h);
        h1  = b(xb);
        err = norm(h1 - h, inf);
        h   = h1
    end % while ()

    h_hat = fft(h);
    h0    = (gam_hat(1) + h_hat(1))/N;
    T     = 1i*kf.*tanh(k*h0);
    T1    = 1i*kf.*coth(k*h0); T1(1) = 0.0;
    T2    = -1i*kf.*csch(k*h0); T2(1) = 0.0;

    gam_x = real(ifft(1i*kf.*k.*gam_hat));
    chi_x = 1 - real(ifft(1i*kf.*k.*(T1.*gam_hat + T2.*h_hat)));
    J     = chi_x.^2 + gam_x.^2;
    J_hat = fft(J);
    J     = real(ifft(kf.*J_hat));
    phi_x = real(ifft(1i*kf.*k.*phi_hat));
    psi_x = real(ifft(1i*k.*T.*phi_hat));

    psixJ = psi_x./J; psixJ_hat = fft(psixJ);
    TpsiJ = real(ifft(T1.*psixJ_hat));
    q     = mean(chi_x.*TpsiJ + gam_x.*psixJ);

    rhs(1:N)     = gam_x.*(TpsiJ + q) - chi_x.*psixJ;
    rhs(N+1:end) = -g*gamma + 0.5*(psi_x.^2 - phi_x.^2)./J + phi_x.*TpsiJ;

end % RHS ()