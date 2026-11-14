% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Frequency of exact steady waves over a flat bed, as the solver computes it. In
% units where g = 1 and the depth is h = 1, Stokes waves of relative depth kh =
% 0.5, 1 and 2 and steepness kA from 0.001 to 0.1 (A the crest elevation above
% the mean level) are computed by tests/steady_wave.m and evolved by RHS.m with
% ode45 over four periods. The frequency is measured from the phase of the
% fundamental Fourier mode of the elevation at a fixed label, sixteen times per
% period, and compared with the exact frequency k c. The script also writes the
% transport constant q of the steady wave (the mean of x_xi T(psi_xi/J) +
% gamma_xi psi_xi/J), from which the frequency error -(3/2) q/c of version 1.0
% is predicted.
%
% The variable codedir may name the code folder of another version and tag the
% suffix of the file written, as in example_data.m:
%     frequency_data                                  (file frequency_v2.csv)
%     codedir = '../../../Euler_bottom-v1.0/code'; tag = 'v1'; frequency_data
% (see example_data.m for the copy of version 1.0).

here = fileparts(mfilename('fullpath'));
if (~exist('codedir', 'var')), codedir = fullfile(here, '..', '..', 'code'); end
if (~exist('tag', 'var')), tag = 'v2'; end
addpath(fullfile(here, '..', '..', 'tests'));
addpath(fullfile(codedir, 'methods'));
addpath(fullfile(codedir, 'problem'));

global a0 g h1 h2 k kf tol xi L L1 Nvec N N2
g  = 1.0; h1 = 1.0; h2 = 1.0; L = 2*pi; L1 = pi;
N  = 256; N2 = N/2; Nvec = 1:N2; xi = (1-N2:N2)'*(2*L/N);
k  = [0:N2 1-N2:-1]'*pi/L;
kf = exp(-5.0*(abs(k)/(0.75*max(k))).^19.0);
tol = 1e-12; a0 = 0;
opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-12);

out = zeros(0, 8);
for k0 = [0.5 1 2]
    for kA = [0.001 0.002 0.005 0.01 0.02 0.05 0.1]
        M = 32;
        [a, c, d] = steady_wave(g, h1, k0, M, kA/k0);
        nk  = (0:M)'*k0;
        cth = [0; coth(nk(2:end)*d)];
        gam = cos(xi*nk')*a;
        phi = c*sin(xi*nk')*(a.*cth);
        % the transport constant of the steady wave, with the formulas of RHS.m
        T1    = 1i*kf.*coth(k*d); T1(1) = 0;
        gam_x = real(ifft(1i*kf.*k.*fft(gam)));
        chi_x = 1 - real(ifft(1i*kf.*k.*T1.*fft(gam)));
        J     = real(ifft(kf.*fft(chi_x.^2 + gam_x.^2)));
        psi_x = real(ifft(1i*k.*(1i*kf.*tanh(k*d)).*fft(phi)));
        q     = mean(chi_x.*real(ifft(T1.*fft(psi_x./J))) + gam_x.*psi_x./J);
        % evolution over four periods, sixteen outputs per period
        T  = 2*pi/(k0*c); P = 4;
        [t, V] = ode45(@RHS, (0:16*P)*T/16, [gam; phi], opts);
        m  = round(k0*L/pi) + 1;                  % index of the fundamental mode
        G  = fft(V(:, 1:N).');
        th = unwrap(angle(G(m, :)));
        om = -(th(end) - th(1))/t(end);
        out(end+1, :) = [k0*h1, kA, c, q, k0*c, om, om/(k0*c) - 1, -1.5*q/c];
        fprintf('%s: kh = %.1f, kA = %.3f: q = %+.3e, frequency error %+.2e (predicted for version 1.0: %+.2e)\n', ...
                tag, k0*h1, kA, q, om/(k0*c) - 1, -1.5*q/c);
    end
end

fid = fopen(fullfile(here, ['frequency_' tag '.csv']), 'w');
fprintf(fid, 'kh,kA,c,q,omega_exact,omega,error,predicted_v1\n');
fprintf(fid, '%.2f,%.4f,%.15f,%.9e,%.15f,%.15f,%.6e,%.6e\n', out.');
fclose(fid);
rmpath(fullfile(codedir, 'methods'));
rmpath(fullfile(codedir, 'problem'));
clear codedir tag
