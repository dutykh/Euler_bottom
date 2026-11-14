% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Steady periodic gravity wave (Stokes wave) over a flat bed of depth h, in the
% conformal variables of the solver. The elevation of the free surface at the
% label xi is gamma(xi) = sum_{n=0}^{M} a(n+1) cos(n k0 xi), with its crest at
% xi = 0. In the frame moving with the wave at the phase speed c the free
% surface is a streamline, on which the Bernoulli condition reads
%
%     c^2/(2 J) + g gamma = B,    J = x_xi^2 + gamma_xi^2,    x_xi = 1 + K gamma,
%
% where K is the Fourier multiplier |k| coth(|k| d) and d = h + a(1) the width of
% the conformal strip. The mean level of the free surface is placed at y = 0
% (mean of gamma x_xi equal to zero), and the crest elevation is gamma(0) = A.
% The unknowns a, c^2 and B solve these M+3 equations, the first at M+1 points
% of half a wavelength, by Newton's method from the linear wave, continued in
% the amplitude over eight steps.
%
% [a, c, d, res] = steady_wave(g, h, k0, M, A) returns the coefficients a, the
% phase speed c, the strip width d and the largest residual res.

function [a, c, d, res] = steady_wave (g, h, k0, M, A)

    n  = (0:M)';
    xj = (0:M)'*pi/(k0*M);                    % crest to trough
    C  = cos(xj*(n'*k0));
    S  = sin(xj*(n'*k0));

    c2 = g/k0*tanh(k0*h);                     % the linear wave
    X  = [zeros(M+1, 1); c2; c2/2];           % unknowns: a(1..M+1), c^2, B
    for step = 1:8
        As = A*step/8;
        if (step == 1)
            X(2) = As;
        end
        for it = 1:50
            R = residual(X, As, g, h, k0, n, C, S);
            if (max(abs(R)) < 1e-14)
                break
            end
            Jm = zeros(M+3);
            for j = 1:M+3
                dx = 1e-7*max(1, abs(X(j)));
                Xp = X; Xp(j) = Xp(j) + dx;
                Xm = X; Xm(j) = Xm(j) - dx;
                Jm(:, j) = (residual(Xp, As, g, h, k0, n, C, S) - residual(Xm, As, g, h, k0, n, C, S))/(2*dx);
            end
            X = X - Jm\R;
        end
    end

    a   = X(1:M+1);
    c   = sqrt(X(M+2));
    d   = h + a(1);
    res = max(abs(residual(X, A, g, h, k0, n, C, S)));

end % steady_wave ()

function R = residual (X, A, g, h, k0, n, C, S)

    M  = numel(n) - 1;
    a  = X(1:M+1);
    c2 = X(M+2);
    B  = X(M+3);
    d  = h + a(1);
    K  = [0; n(2:end)*k0.*coth(n(2:end)*k0*d)];  % multiplier of x_xi - 1
    gam   = C*a;
    gam_x = -S*(n*k0.*a);
    chi_x = 1 + C*(K.*a);
    J     = chi_x.^2 + gam_x.^2;
    % mean of gamma x_xi over a wavelength, by orthogonality of the cosines:
    mlev  = a(1) + 0.5*sum(a(2:end).^2.*K(2:end));
    R     = [c2./(2*J) + g*gam - B; mlev; sum(a) - A];

end % residual ()
