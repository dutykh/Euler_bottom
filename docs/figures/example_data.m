% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% The example of runEuler.m (a wave packet over a depth step), computed as
% runEuler.m computes it (same parameters, initial condition, integrator and
% tolerances, outputs every 10 s up to 200 s), without its plots. For each
% output time the script writes the packet centroid, the mass and the energy;
% for t = 0, 100 and 200 s the free surface, and at t = 0 the bed, in the
% physical plane. MATLAB integrates with ode113, as runEuler.m; GNU Octave,
% which has no ode113, with ode45.
%
% The variable codedir may name the code folder of another version (by default
% the one of this repository) and tag the suffix of the files written:
%     example_data                                   (files example_v2_*.csv)
%     codedir = '../../../Euler_bottom-v1.0/code'; tag = 'v1'; example_data
% Version 1.0 is obtained beside this one by: git worktree add ../Euler_bottom-v1.0 v1.0,
% with a semicolon added at the end of line 33 of its code/methods/RHS.m, without which it
% prints the bed at every evaluation.
% For version 1.0 the bed written is the one its conformal map describes,
% y = h - 2 <h> (its formulas take the depth h for the elevation of the bed).

here = fileparts(mfilename('fullpath'));
if (~exist('codedir', 'var')), codedir = fullfile(here, '..', '..', 'code'); end
if (~exist('tag', 'var')), tag = 'v2'; end
addpath(fullfile(codedir, 'methods'));
addpath(fullfile(codedir, 'problem'));

global a0 g h1 h2 k kf tol xi L L1 Nvec N N2
g    = 9.8;  L1 = 5000.0; L2 = 5000.0; L = L1 + L2;
N    = 4096; N2 = N/2; dxi = 2*L/N; xi = (1-N2:N2)'*dxi; Nvec = 1:N2; tol = 1e-9;
k    = [0:N2 1-N2:-1]'*pi/L; kmax = max(k);
kf   = exp(-5.0*(abs(k)/(0.75*kmax)).^19.0);
h1   = 50.0; h2 = 25.0; a0 = 1.0;
v    = InitCond();
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
if (exist('ode113')), solver = @ode113; else, solver = @ode45; end

t = 0.0; Tf = 200; dtw = 10.0;
out = zeros(0, 5); profiles = zeros(0, 3);
tic;
while (true)
    % the converged bed and the abscissas, as RHS.m and Plot.m compute them
    gam = v(1:N); phi = v(N+1:end);
    h = b(xi); err = inf;
    while (err > tol)
        xb = x_bot(gam, h); hn = b(xb); err = norm(hn - h, inf); h = hn;
    end
    xb = x_bot(gam, h); xs = x(gam, h);
    x_xi = 1 + real(ifft(1i*k.*fft(xs - xi)));          % dx/dxi
    d = mean(gam) + mean(h);                             % width of the strip
    G = k.*tanh(k*d);                                    % Dirichlet-Neumann operator
    M = sum(gam.*x_xi)*dxi;                              % mass above y = 0
    E = 0.5*sum(phi.*real(ifft(G.*fft(phi))))*dxi + 0.5*g*sum(gam.^2.*x_xi)*dxi;
    w = gam(Nvec).^2.*x_xi(Nvec);
    Xc = sum((xs(Nvec) + L).*w)/sum(w);                  % centroid of the packet of the first half
    asym = max(abs(v - [flipud(v(1:N)); flipud(v(N+1:end))]));
    out(end+1, :) = [t, Xc, M, E, asym];
    if (any(abs(t - [0 100 200]) < 1e-9))
        profiles = [profiles; t*ones(N, 1), xs + L, gam];
    end
    if (t == 0)
        if (strcmp(tag, 'v1')), ybed = h - 2*mean(h); else, ybed = -h; end
        bed = [xb + L, ybed];
    end
    fprintf('%s: t = %5.1f s, centroid %.3f m, mass %+.6e, energy %.9e (%.0f s)\n', tag, t, Xc, M, E, toc);
    if (t >= Tf), break; end
    [~, V] = solver(@RHS, [0, dtw], v, options);
    v = V(end, :).';
    t = t + dtw;
end

fid = fopen(fullfile(here, ['example_' tag '_series.csv']), 'w');
fprintf(fid, 't,centroid,mass,energy,asymmetry\n');
fprintf(fid, '%.1f,%.9f,%.15e,%.15e,%.3e\n', out.');
fclose(fid);
fid = fopen(fullfile(here, ['example_' tag '_profiles.csv']), 'w');
fprintf(fid, 't,x,eta\n');
fprintf(fid, '%.1f,%.9f,%.12e\n', profiles.');
fclose(fid);
fid = fopen(fullfile(here, ['example_' tag '_bed.csv']), 'w');
fprintf(fid, 'x,y\n');
fprintf(fid, '%.9f,%.12e\n', bed.');
fclose(fid);
rmpath(fullfile(codedir, 'methods'));
rmpath(fullfile(codedir, 'problem'));
clear codedir tag
