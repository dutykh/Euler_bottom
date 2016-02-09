%%% -------------------------------------------------- %%%
%%% Plotting function to display the solution          %%%
%%% -------------------------------------------------- %%%
%%% Last modified: 09/02/2016                          %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function Plot (v, t)

    % Global variables:
    global a0 g kf xi L Nvec N N2

    % retrieve the free surface elevation
    eta0 = v(Nvec);
    phi0 = v(N + Nvec);

    xf = x(v(1:N), b(xi));
    xf = xf(Nvec);

    xb = x_bot(v(1:N), b(xi));
    h0 = b(x_bot(v(1:N), b(xb)));
    xb = xb(Nvec);

    % We plot the free surface elevation:
    subplot(4,1,1)
    p1 = plot(xf + L, eta0, 'LineWidth', 1.2); grid off, hold on
    hb = b(xb); hb = hb(Nvec);
    p2 = plot(xb + L, -hb, 'k-', 'LineWidth', 1.2); hold on

    opts = {'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75]};
    fill_between(xb + L, min(-h0)*ones(size(xb)), -hb, 1, opts{:}); hold on

    opts = {'EdgeColor', 'none', 'FaceColor', [64 179 255]/255, 'FaceAlpha', 0.5};
    hf   = b(xf);
    fill_between(xf + L, -hf(Nvec), eta0, 1, opts{:}); hold off

    xlim([0.0; L]); ylim([min(-h0)-1.0; 1.0+2.0*a0]);
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', 14);
    ylabel('$\eta(x,t)$', 'interpreter', 'latex', 'fontsize', 14);
    legend([p1, p2], {'Free surface', 'Bottom profile'},...
        'location', 'SouthEast');
    title(['Numerical solution of the Euler equations at t = ',...
        num2str(t,'%3.2f')], 'interpreter', 'latex', 'fontsize', 12);

    subplot(4,1,2)
    plot(xf + L, eta0, 'LineWidth', 1.2); grid off, hold off
    xlim([0.0; L]); ylim([-2.0*a0 2.0*a0]);
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', 14);
    ylabel('$\eta(x,t)$', 'interpreter', 'latex', 'fontsize', 14);
    legend('Free surface elevation');

    subplot(4,1,3)
    plot(xf + L, phi0, 'Color', [252, 78, 122]/255, 'LineWidth', 1.2); grid off, hold off
    xlim([0.0; L]); ylim([-2.0*g*a0 2.0*g*a0]);
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', 14);
    ylabel('$\varphi(x,t)$', 'interpreter', 'latex', 'fontsize', 14);
    legend('Velocity potential at free surface');

    subplot(4,1,4)
    eta_hat = fft(v(1:N));
    FPS     = abs(eta_hat/N).^2;
    I       = find(FPS(Nvec) > 0.0); % we remove zeros to avoid problems with Log ()
    semilogy(Nvec(I), FPS(I), 'k-', 'LineWidth', 0.2), grid off, hold on
    semilogy(Nvec, kf(Nvec), 'r--', 'LineWidth', 1.0), hold off
    xlabel('$N$', 'interpreter', 'latex', 'fontsize', 14);
    ylabel('$|\hat{\eta}(k,t)|^2$', 'interpreter', 'latex', 'fontsize', 14);
    legend('Fourier power spectrum', 'De-aliasing function', 'location', 'NorthEast');
    xlim([0.0 1.05*N2]); ylim([1e-21 1]);

    set(gcf, 'Color', 'w');
    drawnow

end % Plot ()