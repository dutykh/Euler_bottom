% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% The anchor label of the conformal map, xi(1) + s with s = -L/N, half-way
% between the last and the first labels of the periodic grid. The bed point of
% this label is held at x = -L: x_bot.m and x.m take it as the origin of the
% abscissas, and the transport constant q of RHS.m keeps it still. Under the
% mirror symmetry v(j) = v(N+1-j) that b.m and InitCond.m build with flipud,
% this label is a mirror point of the grid.
%
% A periodic function u sampled on the labels takes at the anchor label the
% value real(sum(fft(u).*w))/N of its trigonometric interpolant, in which the
% Nyquist mode enters with the weight cos(k s) = 0.

function [s, w] = anchor ()

    global k L N
    persistent s0 w0 L0 N0

    % the weights depend only on L and N, through k:
    if (isempty(w0) || L0 ~= L || N0 ~= N)
        s0 = -L/N;
        w0 = exp(1i*k*s0);
        w0(N/2+1) = 0.0;
        L0 = L;
        N0 = N;
    end % if ()

    s = s0;
    w = w0;

end % anchor ()
