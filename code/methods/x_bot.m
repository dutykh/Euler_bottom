%%% -------------------------------------------------- %%%
%%% Recover the horizontal coordinate on the bottom    %%%
%%% -------------------------------------------------- %%%
%%% Last modified: 08/02/2016                          %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function xb = x_bot (eta, h)
    
    global k xi L N
    
    % Take the Fourier transforms:
    eta_hat = fft(eta);
    h_hat   = fft(h);
    
    % mean water depth
    H       = (eta_hat(1) + h_hat(1))/N;
    
    %%% Define the pseudo-diff operators:
    T1      = 1i*csch(k*H);
    T1(1)   = 0.0;
    T2      = -1i*coth(k*H);
    T2(1)   = 0.0;
    
    xb      = xi - real(ifft(T1.*eta_hat + T2.*h_hat));
    %%% Remove the drift and put back into the box:
    xb      = xb - xb(1) - L;
    
end % x_bot ()