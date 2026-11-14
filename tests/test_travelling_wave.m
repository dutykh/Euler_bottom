% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Test of the evolution equations on an exact steady wave over a flat bed. In
% units where g = 1 and the depth is h = 1, the Stokes wave of wavenumber k0 = 1
% (kh = 1) whose crest rises 0.1 above the mean level is computed to round-off
% by steady_wave.m, on a periodic domain holding two wavelengths. The exact
% solution is the initial wave translated at its phase speed c. The test checks
%
%   (1) that RHS.m moves the elevation and the potential along the labels with
%       one common velocity, as a change of parametrization must, and
%   (2) that the free surface computed by RHS.m and ode45 over one period, mapped
%       to the physical plane by x.m, coincides at every quarter period with the
%       initial wave translated by c t.
%
% The deviations are of the order of the tolerances of the integrator (1e-10 or
% less). Run from MATLAB or GNU Octave:  run('tests/test_travelling_wave.m')

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', 'code', 'methods'));
addpath(fullfile(here, '..', 'code', 'problem'));

global a0 g h1 h2 k kf tol xi L L1 Nvec N N2
g    = 1.0;                          % gravity acceleration
h1   = 1.0; h2 = 1.0;                % b.m then returns the constant depth 1
L    = 2*pi; L1 = pi;                % the domain (-L, L] holds two wavelengths
N    = 256; N2 = N/2; Nvec = 1:N2;
xi   = (1-N2:N2)'*(2*L/N);
k    = [0:N2 1-N2:-1]'*pi/L;
kf   = exp(-5.0*(abs(k)/(0.75*max(k))).^19.0);   % filter of runEuler.m
tol  = 1e-12; a0 = 0;

% the steady wave: gamma = sum a_n cos(n k0 xi), phi = c sum a_n coth(n k0 d) sin(n k0 xi)
k0 = 1.0; A = 0.1; M = 32;
[a, c, d, res] = steady_wave(g, h1, k0, M, A);
nk  = (0:M)'*k0;
cth = [0; coth(nk(2:end)*d)];
gam0  = cos(xi*nk')*a;
phi0  = c*sin(xi*nk')*(a.*cth);
gam_x = -sin(xi*nk')*(a.*nk);
phi_x = c*cos(xi*nk')*(a.*cth.*nk);
v0 = [gam0; phi0];
fprintf('Steady wave: kh = %.1f, crest elevation %.2f, phase speed c = %.12f (residual %.1e)\n', k0*h1, A, c, res);

% (1) one common velocity along the labels for the elevation and the potential
r = RHS(0, v0);
r(N+1:end) = r(N+1:end) - mean(r(N+1:end));      % the Bernoulli constant
s = -(r(1:N)'*gam_x)/(gam_x'*gam_x);
dev1 = max(abs([r(1:N) + s*gam_x; r(N+1:end) + s*phi_x]));
fprintf('(1) the right-hand side moves both fields with one velocity: residual %.1e\n', dev1);

% (2) evolution over one period and comparison in the physical plane
T = 2*pi/(k0*c);
opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-12);
[t, V] = ode45(@RHS, [0 T/4 T/2 3*T/4 T], v0, opts);

% initial abscissas x0(xi) = xs + xi + sum a_n coth(n k0 d) sin(n k0 xi), the shift
% xs being the origin that x.m imposes; the exact elevation at the abscissa X is
% gamma0(xi) where x0(xi) = X, found by Newton's method.
hb = b(xi);
x0 = x(gam0, hb);
xs = x0(1) - xi(1) - sin(xi(1)*nk')*(a.*cth);
dev2 = zeros(numel(t)-1, 1);
for it = 2:numel(t)
    gam = V(it, 1:N)';
    X   = x(gam, hb) - c*t(it);          % abscissas carried back to t = 0
    z   = X - xs;                        % first guess of the labels
    for newton = 1:30
        z = z - (xs + z + sin(z*nk')*(a.*cth) - X)./(1 + cos(z*nk')*(a.*cth.*nk));
    end
    dev2(it-1) = max(abs(gam - cos(z*nk')*a));
    fprintf('(2) t = %.4f (%.2f period): largest deviation from the translated wave %.1e\n', ...
            t(it), t(it)/T, dev2(it-1));
end

if (dev1 < 1e-12 && max(dev2) < 1e-8)
    fprintf('test passed\n');
else
    error('test FAILED: the wave is not translated rigidly at its phase speed');
end
