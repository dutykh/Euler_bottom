% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Test of the conformal map on an exact map. The analytic function
%
%     Z(w) = w + i m + sum_n c_n exp(i n w),    w = xi + i eta,    -d < eta < 0,
%
% maps the strip of width d = 0.8 onto a fluid layer of period 2 pi. The line
% eta = 0 goes to the free surface and the line eta = -d to the bed, whose depth
% varies between 0.55 and 0.87. Since |Z' - 1| < 1 on both lines, Re Z' > 0 in
% the whole strip and the map is one-to-one. Given at each label xi the
% elevation gamma = Im Z(xi) of the free surface and the depth h = -Im Z(xi - i d)
% of the bed, as b.m gives it, the functions x.m and x_bot.m must return the
% abscissas Re Z(xi) and Re Z(xi - i d) of the surface and bed points, in the
% frame in which the bed point of the anchor label (anchor.m) lies at x = -L.
% Both errors are of the order of 1e-15.
%
% Run from MATLAB or GNU Octave:  run('tests/test_conformal_map.m')

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'code', 'methods'));

global k xi L N
L  = pi;                                 % half-period of the layer
N  = 128;                                % number of labels
xi = (1-N/2:N/2)'*(2*L/N);               % labels, the conformal abscissas
k  = [0:N/2 1-N/2:-1]'*pi/L;             % wavenumbers

d  = 0.8;                                % width of the conformal strip
m  = 0.05;                               % mean elevation of the free surface
n  = [-3 -2 -1 1 2 3];                   % modes of the map and their amplitudes:
c  = [0.004-0.002i, -0.010+0.006i, 0.030+0.012i, -0.025+0.020i, 0.015-0.008i, 0.006+0.003i];
Z  = @(w) w + 1i*m + exp(1i*w*n)*c.';    % the map, for a column of points w
dZ = @(w) 1 + exp(1i*w*n)*(1i*n.*c).';   % its derivative

zs    = Z(xi);                           % surface points
zb    = Z(xi - 1i*d);                    % bed points
gam   = imag(zs);                        % elevation of the free surface
depth = -imag(zb);                       % depth of the bed
za    = Z(xi(1) - L/N - 1i*d);           % bed point of the anchor label

fprintf('Exact map: |Z''-1| at most %.2f on the surface and %.2f on the bed\n', ...
        max(abs(dZ(xi) - 1)), max(abs(dZ(xi - 1i*d) - 1)));
fprintf('depth of the bed from %.3f to %.3f, strip width %.3f\n', min(depth), max(depth), d);

err_s = max(abs(x(gam, depth) - (real(zs) - real(za) - L)));
err_b = max(abs(x_bot(gam, depth) - (real(zb) - real(za) - L)));
fprintf('largest error of the surface abscissas (x.m):   %.1e\n', err_s);
fprintf('largest error of the bed abscissas (x_bot.m):   %.1e\n', err_b);

if (max(err_s, err_b) < 1e-12)
    fprintf('test passed\n');
else
    error('test FAILED: the map does not reproduce the exact one');
end
