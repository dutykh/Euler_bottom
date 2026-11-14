%%% -------------------------------------------------- %%%
%%% Bathymetry function: we use a step at the bottom   %%%
%%% -------------------------------------------------- %%%
%%% Last modified: 14/11/2026                          %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function h = b(x)

    global h1 h2 L L1 Nvec

    % step' steepness:
    Delta = 0.005;  % parameter controlling transition between h1 and h2

    % and the bathymetry:
    x0 = x(Nvec);
    h0 = h1 + (h2 - h1)*(0.5 + 0.5*tanh(Delta*(L + x0 - L1)));
    % the second half of the labels holds the mirror image of the first:
    h  = [h0; flipud(h0)];
    
end % b()