%% example_single_run.m
%  Minimal driver: one nominal-point GNLSE evolution + GSL diagnostics.
%
%  Runtime: about 20-40 s on a contemporary laptop, MATLAB R2021a+.
%  No toolboxes required.
%
%  Produces:
%     - on-screen summary (Delta S_tot, eta_GSL, k_c drift, photon-number drift)
%     - figure: spectral evolution, S_hor/S_rad/S_tot vs xi, eta_GSL vs xi
%     - figure file: example_single_run.png (saved next to the script)

clear; clc; close all;

% --- 1. Run the solver ----------------------------------------------------
out = gnlse_dimensionless( ...
    'N_sol',   3.5, ...
    'delta3',  0.02, ...
    'xi_max',  12.0, ...
    'verbose', true);

% --- 2. Headline numbers --------------------------------------------------
dS_tot = out.S_tot(end) - out.S_tot(1);
eta_f  = out.eta(end);
drift  = (out.P_norm(end) - out.P_norm(1)) / out.P_norm(1);

frac_neg = mean(diff(out.S_tot) < 0);

fprintf('\n=========================================================\n');
fprintf('  SOTHE nominal-point reproducibility summary\n');
fprintf('=========================================================\n');
fprintf('  Delta S_tot         = %+8.3f nats   (paper range [0.82, 2.42])\n', dS_tot);
fprintf('  eta_GSL(xi_end)     = %8.3f         (paper range [0.52, 0.62])\n', eta_f);
fprintf('  fraction dS/dxi < 0 = %8.1f %%      (paper qualitative ~50%%)\n', 100*frac_neg);
fprintf('  photon-number drift = %+8.2e        (target |drift| < 1e-3)\n', drift);
fprintf('  k_c drift           = %+8.3f         (truncated GNLSE: ~0)\n', ...
        out.k_c_hist(end) - out.k_c_hist(1));
fprintf('=========================================================\n');

% --- 3. Quick visualisation ----------------------------------------------
fig = figure('Color','w','Position',[60 60 1300 380]);

subplot(1,3,1)
k       = out.k;
imagesc(out.xi, k, 10*log10(abs(fftshift(fft(out.psi_hist,[],1),1)).^2 + 1e-15));
axis xy; colormap(jet); caxis([-60 10]); colorbar;
xlabel('\xi'); ylabel('k (1/T_0)');
title('(a) Spectral evolution');
ylim([-30 30]);

subplot(1,3,2)
plot(out.xi, out.S_hor, 'b--', 'LineWidth', 1.4); hold on
plot(out.xi, out.S_rad, 'r-.', 'LineWidth', 1.4);
plot(out.xi, out.S_tot, 'k-',  'LineWidth', 2.0);
grid on; xlabel('\xi'); ylabel('S (nats)');
legend('S_{hor}','S_{rad}','S_{tot}','Location','best');
title('(b) Spectral entropies');

subplot(1,3,3)
plot(out.xi, out.eta, 'k-', 'LineWidth', 2.0); hold on
yline(0.5, 'r--', 'symmetric-walk baseline');
ylim([0 1]); grid on;
xlabel('\xi'); ylabel('\eta_{GSL}');
title('(c) Coarse-grained GSL efficiency');

sgtitle(sprintf('SOTHE nominal point: N=%.1f, \\delta_3=%.3f, \\xi_{max}=%.1f', ...
                out.params.N_sol, out.params.delta3, out.params.xi_max));

print(fig, 'example_single_run.png', '-dpng', '-r200');
fprintf('\n  Figure saved: example_single_run.png\n');
