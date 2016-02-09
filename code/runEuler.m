%%% -------------------------------------------------- %%%
%%% Pseudo-spectral solver for full Euler equations    %%%
%%% Test-case: wave packet interaction with a step     %%%
%%% -------------------------------------------------- %%%
%%% Last modified: 09/02/2016                          %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

close all,
clear all,
format longE

%%% Some libraries we use:
addpath odebar/
addpath methods/
addpath problem/
addpath export_fig/

% enable FFTW library
fftw('planner', 'hybrid');

%%% Declaration of global constants and variables:
global a0 g h1 h2 k kf tol xi L L1 Nvec N N2

%%% Physical problem parameters:
g  = 9.8;       % gravity acceleration
L1 = 5000.0;    % domain length before step
L2 = 5000.0;    % domain length after step
L  = L1 + L2;   % domain half-length

%%% Numerical parameters:
% (notice that you use effectively only N/2 modes,
% since the domain is repeated two times)
% (You should count also the losses for anti-aliasing)
N    = 4096;    % number of Fourier modes
N2   = N/2;     % half number of points
dxi  = 2*L/N;   % distance between two physical points
xi   = (1-N2:N2)'*dxi; % physical space discretization
Nvec = 1:N2;    % needed to plot the Fourier power spectrum
tol  = 1e-9;    % tolerance parameter in fixed point iterations

%%% Definition of various (pseudo-)differential operators:
k    = [0:N2 1-N2:-1]'*pi/L;    % vector of wavenumbers
kmax = max(k);                  % max wavenumber (to show the spectrum)
Nvec = 1:N2;                    % positive frequencies in Matlab representation

%%% Dealising function:
Lambd= 19.0;
kam  = 0.75*kmax;
kf   = exp(-5.0*(abs(k)/kam).^Lambd);

%%% Domain and bathymetry definitions:
x0 = L + xi(1:N2);          % half-x domain
h1 = 50.0;                  % depth before the step
h2 = 25.0;                  % depth after the step

%%% Definition of initial condition:
% Parameters:
a0 = 1.0;       % wave max amplitude
v  = InitCond();% (problem-dependent)

%%% Time integration parameters and options:
t   = 0.0;  % the discrete time variable
Tf  = 200;  % final simulation time
dtw = 10.0; % display results every dtw seconds

% time-stepping parameters for the ODE solver:
options  = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'OutputFcn', @odetpbar);

%%% Plot the initial condition:
FigHandle = figure(1); % Figure window size
set(FigHandle, 'Renderer', 'zbuffer');
set(FigHandle, 'Position', [900, 184, 1020, 990]);
Plot(v, t);

while (t <= Tf) % the loop in time
    [~, V] = ode113(@RHS, [0, dtw], v, options);
    v      = V(end,:).';
    Plot(v, t);
    t      = t + dtw;
end % while ()