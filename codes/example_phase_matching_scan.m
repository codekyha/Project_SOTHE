%% example_phase_matching_scan.m
%  Reproduces the headline phase-matching scaling: k_RR * delta_3 ~ const.
%  Replaces the withdrawn Boltzmann thermal-tail fit of the original
%  submission (see NOTES.md, Section 2).
%
%  Scans delta_3 in [0.06, 0.10] at fixed N = 3.5 and fits
%      k_RR = c_RR / delta_3
%  through the origin. Paper values: c_RR = 1.01, mean k_RR*delta_3 = 1.02,
%  R^2 = 0.93.
%
%  Runtime: about 4-6 minutes (5 GNLSE runs, ~1 min each).

clear; clc; close all;

delta3_list = [0.06, 0.07, 0.08, 0.09, 0.10];
N_sol       = 3.5;
xi_max      = 8.0;          % shorter than the GSL run; the lobe is set early

k_RR = zeros(size(delta3_list));
for i = 1:numel(delta3_list)
    fprintf('--- run %d/%d: delta_3 = %.3f ---\n', i, numel(delta3_list), delta3_list(i));
    out = gnlse_dimensionless('N_sol', N_sol, 'delta3', delta3_list(i), ...
                              'xi_max', xi_max, 'verbose', false);

    % final-snapshot spectrum
    Spec = abs(fftshift(fft(out.psi_hist(:, end)))).^2;
    k_c  = out.k_c_hist(end);
    k_RR(i) = find_rr_lobe(out.k, Spec, k_c, '+', 2.0);
    fprintf('   k_RR = %.3f   (k_RR * delta_3 = %.3f)\n', k_RR(i), k_RR(i)*delta3_list(i));
end

fit = rr_phase_match_fit(delta3_list, k_RR);

fprintf('\n=========================================================\n');
fprintf('  Phase-matching scan summary\n');
fprintf('=========================================================\n');
fprintf('  c_RR             = %.3f   (paper: 1.01)\n', fit.c_RR);
fprintf('  mean k_RR*delta3 = %.3f   (paper: 1.02)\n', mean(fit.product));
fprintf('  R^2              = %.3f   (paper: 0.93)\n', fit.R2);
fprintf('=========================================================\n');

% --- plot ---
fig = figure('Color','w','Position',[80 80 700 400]);
d3_dense = linspace(min(delta3_list)*0.95, max(delta3_list)*1.05, 200);
plot(d3_dense, fit.c_RR ./ d3_dense, 'b-', 'LineWidth', 2.0); hold on
plot(delta3_list, k_RR, 'ro', 'MarkerSize', 9, 'MarkerFaceColor','w', 'LineWidth', 1.4);
grid on; xlabel('\delta_3'); ylabel('k_{RR}  (1/T_0)');
legend(sprintf('k_{RR} = %.3f / \\delta_3', fit.c_RR), 'measured', 'Location','northeast');
title(sprintf('Cherenkov phase-matching scan: c_{RR}=%.3f, R^2=%.2f', fit.c_RR, fit.R2));

print(fig, 'example_phase_matching_scan.png', '-dpng', '-r200');
fprintf('\n  Figure saved: example_phase_matching_scan.png\n');

% --- export numerical data ---
T = table(delta3_list(:), k_RR(:), k_RR(:).*delta3_list(:), ...
          'VariableNames', {'delta3','k_RR','k_RR_times_delta3'});
writetable(T, 'example_phase_matching_scan.csv');
fprintf('  CSV saved:    example_phase_matching_scan.csv\n');
