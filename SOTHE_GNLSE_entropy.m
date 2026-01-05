%% ========================================================================
%  Optical Entropy, Solitonic Horizon & Hawking Temperature Extraction
%  GNLSE with beta2, beta3, Kerr nonlinearity
%  Demonstrates Generalized Second Law + Effective Hawking Temperature
%
%  Author: H. Oguz
%% ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%% 1. Physical Parameters
%% ------------------------------------------------------------------------

c_light = 299792.458;            % Speed of light [nm/ps]
lambda0 = 800;                  % Central wavelength [nm]
omega0  = 2*pi*c_light/lambda0; % rad/ps

beta2 = -15.0;                  % ps^2/km (anomalous)
beta3 = +0.10;                  % ps^3/km
gamma = 0.1;                    % 1/W/km

N_sol = 3;                      % Soliton order (fission regime)
T0    = 0.050;                  % Pulse width [ps]

%% ------------------------------------------------------------------------
%% 2. Numerical Grid
%% ------------------------------------------------------------------------

n_points = 2^14;
T_window = 10.0;                % ps

dt = T_window / n_points;
t  = (-n_points/2:n_points/2-1).' * dt;

dw = 2*pi / T_window;
w  = (-n_points/2:n_points/2-1).' * dw;

%% ------------------------------------------------------------------------
%% 3. Initial Conditions
%% ------------------------------------------------------------------------

LD = T0^2 / abs(beta2) * 1e-3;          % Dispersion length [km]
P0 = (N_sol^2 * abs(beta2)) / (gamma*T0^2);

U = sqrt(P0) * sech(t/T0);              % Initial soliton

D_op = 1i * (beta2/2 * w.^2 + beta3/6 * w.^3);

%% ------------------------------------------------------------------------
%% 4. Spectral Partition (Thermodynamic Split)
%% ------------------------------------------------------------------------

filter_width = 4 / T0;
mask_sol = exp(-(w/filter_width).^8);
mask_rad = 1 - mask_sol;

%% ------------------------------------------------------------------------
%% 5. Split-Step Fourier Method
%% ------------------------------------------------------------------------

L_prop  = 80 * LD;
n_steps = 12000;
dz      = L_prop / n_steps;

half_step_disp = exp(D_op * dz / 2);

n_saves = 300;
save_idx = round(linspace(1,n_steps,n_saves));
z_save = linspace(0,L_prop,n_saves);

U_save = zeros(n_points,n_saves);
S_rad_hist = zeros(1,n_saves);
S_hor_hist = zeros(1,n_saves);

counter = 1;

fprintf('Running SSFM: %.1f L_D, dz = %.2e km\n',L_prop/LD,dz);

for k = 1:n_steps

    U = ifft( fft(U) .* half_step_disp );
    U = U .* exp(1i*gamma*abs(U).^2*dz);
    U = ifft( fft(U) .* half_step_disp );

    if k == save_idx(counter)

        U_save(:,counter) = U;

        Spec = abs(fftshift(fft(U))).^2;
        W = fftshift(w);
        Msol = fftshift(mask_sol);

        % Horizon entropy proxy
        S_hor_hist(counter) = trapz(W, Spec .* Msol);

        % Radiation entropy
        Spec_rad = Spec .* (1 - Msol);
        Prad = trapz(W, Spec_rad);

        if Prad > 1e-18
            p = Spec_rad / Prad;
            valid = p > 1e-15;
            S_rad_hist(counter) = ...
                -trapz(W(valid), p(valid).*log(p(valid)));
        end

        counter = counter + 1;
        if counter > n_saves, break; end
    end
end

fprintf('Simulation complete.\n');

%% ------------------------------------------------------------------------
%% 6. Robust Hawking Temperature Extraction (FIXED)
%% ------------------------------------------------------------------------
idx_fit = round(0.95*n_saves);          % Use very end of simulation
U_fit   = U_save(:,idx_fit);
Spec    = abs(fftshift(fft(U_fit))).^2;
omega   = fftshift(w);

% --- 1. Locate Resonant Radiation (RR) Peak ---
% We look for the highest peak in the dispersive wave region
% (Assuming beta3 > 0, RR is usually at positive detuning relative to zero)
Msol = fftshift(mask_sol);
Spec_rad = Spec .* (1 - Msol); 

% Find max peak outside soliton bandwidth
[max_val, idx_RR] = max(Spec_rad);
omega_RR = omega(idx_RR);

% --- 2. Define Thermal Tail Region ---
% We scan slightly BEYOND the RR peak to find the exponential decay.
% Reduced buffer to 10% to capture the tail immediately after the peak.
omega_buffer = 0.10 * abs(omega_RR); 

if omega_RR > 0
    % RR is on the right (Blue shift) -> Tail is to the right
    fit_mask = (omega > omega_RR + omega_buffer) & ...
               (Spec > 1e-15*max(Spec)) & ...
               (Spec < 0.5*max_val); % Look below the peak
else
    % RR is on the left (Red shift) -> Tail is to the left
    fit_mask = (omega < omega_RR - omega_buffer) & ...
               (Spec > 1e-15*max(Spec)) & ...
               (Spec < 0.5*max_val);
end

omega_fit = omega(fit_mask);
Spec_fit  = Spec(fit_mask);

% --- 3. Perform Fit Safety Check ---
if numel(omega_fit) < 10
    warning('Still not enough points. Skipping T_H extraction.');
    T_H_eff = NaN;
    lnI = []; 
    p = [];
else
    lnI = log(Spec_fit);
    
    % Linear regression: ln(I) = - (hbar/kB*T) * w + C
    p = polyfit(omega_fit, lnI, 1);
    slope = p(1); % Units: [ps]
    
    % Physical constants
    hbar_J = 1.0545718e-34;
    kB_J   = 1.380649e-23;
    
    % Calculate Temperature
    % Slope is negative for decay. T = -hbar / (kB * slope_SI)
    % Convert slope from ps to s: slope_SI = slope * 1e-12
    slope_SI = slope * 1e-12; 
    
    T_H_eff = abs(hbar_J / (kB_J * slope_SI)); 
end

if ~isnan(T_H_eff)
    fprintf('RR Peak found at: %.2f THz\n', omega_RR/(2*pi));
    fprintf('Effective Hawking Temperature: T_H = %.1f K\n', T_H_eff);
else
    fprintf('Could not extract Temperature (Pulse might not have fissioned yet).\n');
end

%% ------------------------------------------------------------------------
%% 7. Visualization (Safe Plotting)
%% ------------------------------------------------------------------------
S_hor = S_hor_hist / max(S_hor_hist);
S_rad = S_rad_hist / max(S_rad_hist);
S_tot = S_hor + S_rad;

freq_THz = fftshift(w)/(2*pi);
Spectrogram_dB = 10*log10(abs(fftshift(fft(U_save,[],1),1)).^2 + 1e-16);

figure('Color','w','Position',[100 100 1200 400])

% --- (a) Spectral Map ---
subplot(1,3,1)
imagesc(z_save, freq_THz, Spectrogram_dB)
axis xy
ylim([-100 100]) % Adjust zoom if RR is outside
caxis([-40 10])  % Better contrast
xlabel('z (km)')
ylabel('\nu (THz)')
title('(a) Spectral Evolution')
colormap jet
colorbar

% --- (b) Entropy Law ---
subplot(1,3,2)
plot(z_save, S_hor, 'r--', 'LineWidth', 2); hold on
plot(z_save, S_rad, 'b', 'LineWidth', 2)
plot(z_save, S_tot, 'k', 'LineWidth', 1.5)
grid on
xlabel('z (km)')
ylabel('Normalized Entropy')
title('(b) Generalized Second Law')
legend('Horizon S_{BH}', 'Radiation S_{rad}', 'Total S_{tot}', 'Location','best')

% --- (c) Hawking Temperature Fit ---
subplot(1,3,3)
if ~isnan(T_H_eff) && ~isempty(lnI)
    plot(omega_fit/(2*pi), lnI, 'b.', 'MarkerSize', 6); hold on
    plot(omega_fit/(2*pi), polyval(p, omega_fit), 'r', 'LineWidth', 2)
    xlabel('\nu (THz)')
    ylabel('ln(PSD)')
    title(['(c) T_H Fit \approx ' num2str(T_H_eff,'%.0f') ' K'])
    grid on
    legend('RR Tail Data', 'Thermal Fit')
else
    text(0.5, 0.5, 'Fit Failed / No Fission', 'HorizontalAlignment', 'center')
    title('(c) Thermal Tail Analysis')
end

print('Optical_Hawking_Entropy_Fixed','-dpng','-r300')
fprintf('Figure saved.\n');
