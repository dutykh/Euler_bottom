%%% -------------------------------------------------- %%%
%%% Recover the horizontal coorinate on the free surf  %%%
%%% -------------------------------------------------- %%%
%%% Last modified: 14/11/2026                          %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function xf = x (eta, h)

    global k xi L N

    % Take the Fourier transforms:
    eta_hat = fft(eta);
    h_hat   = fft(h);

    % mean water depth
    H       = (eta_hat(1) + h_hat(1))/N;

    %%% Define the pseudo-diff operators:
    T1      = 1i*coth(k*H);
    T1(1)   = 0.0;
    T2      = -1i*csch(k*H);
    T2(1)   = 0.0;

    % h is the depth, so the bed elevation entering the map is -h:
    xf      = xi - real(ifft(T1.*eta_hat - T2.*h_hat));

    %%% The abscissa of the bed point of the anchor label (see anchor.m),
    %%% with the operators of x_bot:
    B1      = 1i*csch(k*H);
    B1(1)   = 0.0;
    B2      = -1i*coth(k*H);
    B2(1)   = 0.0;
    [s, w]  = anchor();
    xs      = xi(1) + s - real(sum((B1.*eta_hat - B2.*h_hat).*w))/N;

    %%% Remove the drift and put back into the box with the origin of
    %%% x_bot, which holds the bed point of the anchor label at x = -L:
    xf      = xf - xs - L;

end % x ()